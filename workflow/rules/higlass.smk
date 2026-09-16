################################################################################
# higlass.smk — Prepare a HiGlass-ingest manifest.
#
# HiGlass (Kerpedjiev et al. 2018) is a web-based interactive Hi-C browser.  The
# `.mcool` files cooler produces are HiGlass-native — no conversion needed.
#
# This rule does two small things:
#   1. Symlinks all per-sample .mcool files into one flat directory
#      (results/higlass/data/) so a local HiGlass server can bind-mount it as
#      its data root in one shot.
#   2. Writes a YAML "manifest" describing the tracks (uid, filetype, datatype)
#      so the launcher script can auto-ingest them.
#
# Spinning up the HiGlass server itself is a host-side concern; see
#   scripts/launch_higlass.sh
# for a one-liner that runs the official higlass/higlass-docker image, mounts
# the data dir, and ingests every .mcool in it.
#
# Gated by config.higlass.enabled: true.
#
# RSchwarz 2026-06-10
################################################################################

rule prep_higlass:
    """
    Build the HiGlass data directory (symlinks of .mcool files) + manifest YAML.
    """
    input:
        mcools = expand(str(RESULTS / "matrix/cool/{sample}.mcool"),
                        sample=SAMPLES),
        cs     = CHROMSIZES,
    output:
        manifest = str(RESULTS / "higlass/manifest.yaml"),
        readme   = str(RESULTS / "higlass/README.md"),
    params:
        datadir = str(RESULTS / "higlass/data"),
        port    = config.get("higlass", {}).get("port", 8989),
        launcher = str(SCRIPTS / "launch_higlass.sh"),
    threads: 1
    resources:
        mem_mb = 1000, runtime = 10
    run:
        import os, textwrap, datetime
        os.makedirs(params.datadir, exist_ok=True)
        # Symlink each .mcool into the data directory
        manifest_entries = []
        for m in input.mcools:
            sample = os.path.basename(m).replace(".mcool", "")
            link = os.path.join(params.datadir, f"{sample}.mcool")
            try:
                if os.path.islink(link) or os.path.exists(link):
                    os.remove(link)
                os.symlink(os.path.realpath(m), link)
            except OSError as e:
                print(f"[prep_higlass] warning: could not symlink {sample}: {e}")
            manifest_entries.append({
                "uid":         sample,
                "filename":    f"{sample}.mcool",
                "filetype":    "cooler",
                "datatype":    "matrix",
            })
        # Also symlink the chrom.sizes for the HiGlass chromosomes track
        cs_link = os.path.join(params.datadir, "chrom_sizes.tsv")
        try:
            if os.path.islink(cs_link) or os.path.exists(cs_link):
                os.remove(cs_link)
            os.symlink(os.path.realpath(input.cs), cs_link)
        except OSError:
            pass
        # Write the manifest YAML
        with open(output.manifest, "w") as fh:
            fh.write("# HiGlass ingest manifest — auto-generated\n")
            fh.write(f"# generated: {datetime.date.today()}\n")
            fh.write("tracks:\n")
            for e in manifest_entries:
                fh.write(f"  - uid: {e['uid']}\n")
                fh.write(f"    filename: {e['filename']}\n")
                fh.write(f"    filetype: {e['filetype']}\n")
                fh.write(f"    datatype: {e['datatype']}\n")
            fh.write("chromsizes:\n")
            fh.write(f"  filename: chrom_sizes.tsv\n")
            fh.write(f"  uid: chromsizes-{config['genome']}\n")
        # Write the README
        with open(output.readme, "w") as fh:
            fh.write(textwrap.dedent(f"""\
                # HiGlass viewing for this Hi-C run

                The `data/` subdirectory contains symlinks to every `{config['genome']}` `.mcool`
                produced by the pipeline. The `manifest.yaml` lists each track ready
                for HiGlass ingest.

                ## To launch a local HiGlass server

                The helper script does it in one line. Requires Docker or Podman:

                ```bash
                bash {params.launcher} \\
                    --data {params.datadir} \\
                    --manifest {output.manifest} \\
                    --port {params.port}
                ```

                It will:
                  1. `docker pull higlass/higlass-docker:latest` (one-time, ~1 GB).
                  2. Start the container with the data dir bind-mounted.
                  3. Ingest every track in the manifest.
                  4. Print the local URL (default http://localhost:{params.port}).

                ## Or — point an existing HiGlass server at the manifest

                If you already have a HiGlass instance, you can ingest from this
                manifest directly:

                ```bash
                python /path/to/higlass-server/manage.py ingest_tileset \\
                    --filename data/<sample>.mcool \\
                    --filetype cooler --datatype matrix --uid <sample>
                ```

                ## Files in this directory

                - `manifest.yaml`  — list of every track ready for ingest
                - `data/*.mcool`   — symlinks to the per-sample matrices
                - `data/chrom_sizes.tsv` — symlink to the genome's chrom.sizes
                """))
