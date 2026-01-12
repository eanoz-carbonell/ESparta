# -*- coding: utf-8 -*-
__author__ = 'benkjohnson & ErnestoAC & ChatGPT' 

import os
import subprocess
import glob
from shutil import copy
import pysam
import check_dependencies_mac
import qc_analysis

class Mapping_and_Counting(object):
    def __init__(self):
        return

    def fix_strand_bam(self, bam_in, bam_out):
        """
        Invert the strand of all reads in a BAM file for reverse stranded libraries.
        Compatible con Python 2.7
        """
        print("Fixing strand for reverse stranded library: {}".format(bam_out))
        infile = pysam.AlignmentFile(bam_in, "rb")
        outfile = pysam.AlignmentFile(bam_out, "wb", template=infile)

        for read in infile:
            read.is_reverse = not read.is_reverse
            outfile.write(read)

        infile.close()
        outfile.close()

    def bowtie(self, datalocation, analysislocation, options):
        """
        Run Bowtie2 for single-end reads.
        Compatible con Python 2.7
        """
        cd = check_dependencies_mac.CheckDependencies()
        qc = qc_analysis.QC_analysis()
        gff, genref = qc.findreferencefiles(datalocation)

        # Copy reference files
        copy(genref, os.path.join(analysislocation, 'Bowtie'))
        genrefname = os.path.basename(genref)
        copy(gff, os.path.join(analysislocation, 'HTSeq'))

        # Prepare Bowtie2 directory
        os.chdir(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting"))
        if not os.path.lexists(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting", "bowtie2")):
            raise Exception("Bowtie2 binary folder not found. Install Bowtie2 or place it inside Mapping_and_counting/bowtie2")
        os.chdir(os.path.join(cd.getpwd(), "bowtie2"))

        # Decompress/copy FASTQ files
        for file in os.listdir(os.path.join(analysislocation, "QC")):
            ext = file.split(".")[-1]
            infile = os.path.join(analysislocation, "QC", file)
            outfile = os.path.join(analysislocation, "Bowtie", os.path.splitext(file)[0])
            if ext == "gz":
                subprocess.Popen("gunzip -c {0} > {1}".format(infile, outfile), shell=True).wait()
            else:
                copy(infile, os.path.join(analysislocation, "Bowtie"))

        # Cleanup QC folder
        if options.cleanup:
            for file in os.listdir(os.path.join(analysislocation, "QC")):
                ext = file.split(".")[-1]
                if ext in ["gz", "fq", "fastq"]:
                    subprocess.Popen("rm {0}".format(os.path.join(analysislocation, 'QC', file)), shell=True).wait()

        # -----------  INDEX BOWTIE2  --------------
        print("Building Bowtie2 index from the reference genome")
        refpath = os.path.join(analysislocation, "Bowtie", genrefname)
        refprefix = os.path.splitext(refpath)[0]

        cmd_index = "./bowtie2-build {0} {1}".format(refpath, refprefix)
        if not options.verbose:
            cmd_index += " > /dev/null 2>&1"
        subprocess.Popen(cmd_index, shell=True).wait()

        # Move index files
        for idx in glob.glob("*.bt2*"):
            copy(idx, os.path.join(analysislocation, "Bowtie"))

        # -----------   MAPPING BOWTIE2 ------------
        print("Mapping reads with Bowtie2")
        ref = os.path.join(analysislocation, "Bowtie", os.path.splitext(genrefname)[0])

        for file in os.listdir(os.path.join(analysislocation, "Bowtie")):
            ext = os.path.splitext(file)[1]
            if ext in [".fq", ".fastq"]:
                fname = os.path.splitext(file)[0]
                strippedfile = fname[len('trimmed'):] if fname.startswith("trimmed") else fname
                infile = os.path.join(analysislocation, "Bowtie", file)
                outfile = os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam")

                cmd = "./bowtie2 -x {0} -U {1} -S {2} --threads {3}".format(ref, infile, outfile, options.threads)
                if options.mismatch:
                    cmd += " -N {0}".format(options.mismatch)
                if options.otherbowtieoptions:
                    cmd += " {0}".format(options.otherbowtieoptions)
                if not options.verbose:
                    cmd += " --quiet"

                print("Running Bowtie2:", cmd)
                subprocess.Popen(cmd, shell=True).wait()

        return

    def htseq(self, analysislocation, options):
        """
        featureCounts + strand fix for single-end reverse stranded libraries.
        Compatible con Python 2.7
        """
        cd = check_dependencies_mac.CheckDependencies()

        # Buscar archivo GTF
        gff_files = glob.glob(os.path.join(analysislocation, "HTSeq") + "/*.gtf") + \
                    glob.glob(os.path.join(analysislocation, "HTSeq") + "/*.gff*")
        if not gff_files:
            raise Exception("No GFF/GTF file found in HTSeq folder")
        gtf = gff_files[0]

        cont = raw_input("Vas a ejecutar featureCounts. ¿Deseas continuar? (y/n): ")
        if cont.lower() not in ["y", "yes"]:
            print("Ejecución detenida por el usuario.")
            return

        print("Counting gene features with featureCounts")

        fc_bin = "/Users/Ernesto/subread-2.1.1-macOS-x86_64/bin/featureCounts"

        for mapfile in os.listdir(os.path.join(analysislocation, "Bowtie")):
            ext = os.path.splitext(mapfile)[1]
            if ext == ".sam":
                samfile = os.path.join(analysislocation, "Bowtie", mapfile)

                # Ordenar → BAM
                sorted_bam = samfile.replace(".sam", ".sorted.bam")
                sort_cmd = "samtools sort -o {0} {1}".format(sorted_bam, samfile)
                print("Sorting for featureCounts:", sort_cmd)
                subprocess.Popen(sort_cmd, shell=True).wait()

                # Index
                subprocess.Popen("samtools index {0}".format(sorted_bam), shell=True).wait()

                # Fix strand for reverse stranded library
                stranded_bam = sorted_bam.replace(".bam", ".stranded.bam")
                self.fix_strand_bam(sorted_bam, stranded_bam)
                subprocess.Popen("samtools index {0}".format(stranded_bam), shell=True).wait()

                # featureCounts
                fname = os.path.splitext(mapfile)[0]
                stripped = fname[len("align"):]
                fc_output = os.path.join(analysislocation, "HTSeq", "fc" + stripped + ".txt")

                cmd = (
                    "{0} -a {1} -o {2} -t exon -g gene_id -Q 0 -s 2 {3}".format(
                        fc_bin, gtf, fc_output, stranded_bam
                    )
                )
                print("Running featureCounts:", cmd)
                subprocess.Popen(cmd, shell=True).wait()

                # Convert featureCounts → .sam estilo SPARTA
                outfile = os.path.join(analysislocation, "HTSeq", "map" + stripped + ".sam")
                convert_cmd = "awk 'BEGIN{{FS=\"\\t\"}} NR>2 {{print $1\"\\t\"$7}}' {0} > {1}".format(fc_output, outfile)
                subprocess.Popen(convert_cmd, shell=True).wait()

        # Cleanup
        if options.cleanup:
            for file in os.listdir(os.path.join(analysislocation, "Bowtie")):
                subprocess.Popen("rm {0}".format(os.path.join(analysislocation, 'Bowtie', file)), shell=True).wait()

        return
