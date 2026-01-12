__author__ = 'benkjohnson & ErnestoAC & ChatGPT'

import os
import subprocess
import numpy as np

class DifferentialExpression(object):
    def __init__(self):
        return

    def removenoncountdata(self, analysislocation):
        """Remove lines starting with '__' from HTSeq .sam count files."""
        countspath = os.path.join(analysislocation, 'HTSeq')
        depath = os.path.join(analysislocation, 'DEanalysis')
        countdata = os.listdir(countspath)
        for countfile in countdata:
            extension = os.path.splitext(countfile)[1]
            if extension == ".sam":
                with open(os.path.join(countspath, countfile), "r") as counts:
                    lines = counts.readlines()
                with open(os.path.join(depath, countfile), "w") as newcounts:
                    for line in lines:
                        if not line.startswith("__"):
                            newcounts.write(line)
        return

    def getuserinput(self, analysislocation):
        """Take user input for specifying the conditions to be tested."""
        depath = os.path.join(analysislocation, 'DEanalysis')
        print "SPARTA has these files:"
        number = 1
        for f in os.listdir(depath):
            print "{0}) {1}".format(number, f)
            number += 1

        moveforward = False
        while not moveforward:
            conditionnumber = raw_input("How many conditions are there?: ")
            if conditionnumber in ['', '0', '1']:
                print "You can't compare a condition to itself."
            else:
                try:
                    int(conditionnumber)
                    usersure = raw_input("Are you sure that's how many conditions you would like to compare? (y/n): ")
                    if usersure.upper() == "Y":
                        moveforward = True
                except Exception:
                    print "Please enter a number (i.e. 2)"

        conditionnumber = int(conditionnumber)
        conditioncounter = 1
        with open(os.path.join(analysislocation, 'DEanalysis', 'conditions_input.txt'), "w") as conditions_input:
            while conditioncounter != conditionnumber+1:
                if conditioncounter == 1:
                    conditions_input.write("Reference_Condition_Files:\n")
                else:
                    conditions_input.write("Experimental_Condition_{0}_Files:\n".format(conditioncounter))
                conditioncounter += 1

        # example file
        exampleconditionnumber = 2
        exampleconditioncounter = 1
        with open(os.path.join(analysislocation, 'DEanalysis', 'conditions_input_example.txt'), "w") as exampleconditions_input:
            while exampleconditioncounter != exampleconditionnumber+1:
                if exampleconditioncounter == 1:
                    exampleconditions_input.write("Reference_Condition_Files: mapWT_A.sam, mapWT_B.sam\n")
                else:
                    exampleconditions_input.write("Experimental_Condition_{0}_Files: mapmutant2466_A.sam, mapmutant2466_B.sam\n".format(exampleconditioncounter))
                exampleconditioncounter += 1

        print "Now edit 'conditions_input.txt' to specify which files belong to each condition, with replicates separated by commas."

        moveforward = False
        conditions_list = []
        while not moveforward:
            proceedanswer = raw_input("Once you have entered the file names, hit Enter: ")
            if proceedanswer == '':
                with open(os.path.join(analysislocation, 'DEanalysis', 'conditions_input.txt'), "r") as user_condition_input:
                    lines = user_condition_input.readlines()
                    try:
                        for line in lines:
                            if line.strip() != '':
                                parts = line.strip().split(':')[1].split(',')
                                conditions_list.append([f.strip() for f in parts])
                        moveforward = True
                    except Exception:
                        print "Error parsing file. Make sure replicates are separated by commas."
        return conditions_list

    def generatecontrasts(self, contrastrow):
        if contrastrow <= 2:
            return None
        contrastcounter = 2
        contrastcol = 0
        while contrastcounter != contrastrow:
            contrastcol += contrastrow - contrastcounter
            contrastcounter += 1
        contrast = np.zeros((contrastcol, contrastrow)).astype(int)
        contrastcounter = 2
        rowindex = 0
        while contrastcounter != contrastrow:
            iterationnumber = contrastrow - contrastcounter
            for itercounter in range(iterationnumber):
                contrast[rowindex][contrastcounter-1] = 1
                contrast[rowindex][(contrastcounter-1)+(itercounter+1)] = -1
                rowindex += 1
            contrastcounter += 1
        return contrast

    def writeRscript(self, analysislocation, conditions_list):
        """Write an R script compatible with modern edgeR and R 4.5+."""
        contrastrow = len(conditions_list)
        contrast = self.generatecontrasts(contrastrow)

        rscript_path = os.path.join(analysislocation, 'DEanalysis', 'DEexpression.r')
        with open(rscript_path, "w") as de_expression:
            de_expression.write("if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n")
            de_expression.write("BiocManager::install('edgeR', ask=FALSE, update=FALSE)\n")
            de_expression.write("library(edgeR)\nlibrary(limma)\n")

            Rvarlst = []
            for t_index, treatments in enumerate(conditions_list):
                repnumcounter = 1
                for varnames in treatments:
                    if varnames != '':
                        varnames = varnames.strip()
                        varname = os.path.splitext(varnames)[0]
                        varname = varname[len('map'):]
                        filepathname = os.path.join(analysislocation, 'DEanalysis', varnames)
                        de_expression.write("{0}_{1} <- read.table('{2}', row.names=1)\n".format(varname, repnumcounter, filepathname))
                        de_expression.write("colnames({0}_{1}) <- '{0}_{1}'\n".format(varname, repnumcounter))
                        Rvarlst.append("{0}_{1}".format(varname, repnumcounter))
                        repnumcounter += 1

            de_expression.write("alldata <- cbind({0})\n".format(",".join(Rvarlst)))

            Rgrouplst = []
            repnum = 1
            for i, treatments in enumerate(conditions_list):
                for group in treatments:
                    if group != '':
                        Rgrouplst.append(str(i+1))
            de_expression.write("group <- factor(c({0}))\n".format(",".join(Rgrouplst)))
            de_expression.write("y <- DGEList(counts=alldata, group=group)\n")
            halfofsamples = len(Rvarlst)/2
            de_expression.write("keep <- rowSums(cpm(y)>2)>={0}\n".format(halfofsamples))
            de_expression.write("y <- y[keep,]\n")
            de_expression.write("y$samples$lib.size <- colSums(y$counts)\n")
            de_expression.write("y <- calcNormFactors(y)\n")
            de_expression.write("png('{0}')\n".format(os.path.join(analysislocation, 'DEanalysis', 'MDSplot.png')))
            de_expression.write("mymdsobj <- plotMDS(y)\n")
            de_expression.write("dev.off()\n")
            de_expression.write("design <- model.matrix(~group)\n")
            de_expression.write("y <- estimateGLMCommonDisp(y, design, verbose=TRUE)\n")
            de_expression.write("y <- estimateGLMTrendedDisp(y, design)\n")
            de_expression.write("y <- estimateGLMTagwiseDisp(y, design)\n")
            de_expression.write("png('{0}')\n".format(os.path.join(analysislocation, 'DEanalysis', 'BCVplot.png')))
            de_expression.write("plotBCV(y)\n")
            de_expression.write("dev.off()\n")
            de_expression.write("fit <- glmFit(y, design)\n")

            # Comparisons to reference
            for condval in range(2, len(conditions_list)+1):
                de_expression.write("lrt <- glmLRT(fit, coef={0})\n".format(condval))
                de_expression.write("FDR <- p.adjust(lrt$table$PValue, method='BH')\n")
                de_expression.write("summary(dt <- decideTests(lrt))\n")
                de_expression.write("isDE <- as.logical(dt)\n")
                de_expression.write("DEnames <- rownames(y)[isDE]\n")
                de_expression.write("png('{0}')\n".format(os.path.join(analysislocation, 'DEanalysis', 'ReferenceCondvsExpCond{0}_scatter.png'.format(condval))))
                de_expression.write("plotSmear(lrt, de.tags=DEnames, cex=0.5)\n")
                de_expression.write("abline(h=c(-1,1), col='blue')\n")
                de_expression.write("dev.off()\n")
                de_expression.write("outfile <- cbind(cpm(y), lrt$table, FDR)\n")
                de_expression.write("write.csv(outfile, file='{0}')\n".format(os.path.join(analysislocation, 'DEanalysis', 'ReferenceCondvsExpCond{0}_DE.csv'.format(condval))))

            # Contrasts between experimental conditions
            if contrast is not None:
                for contr in contrast:
                    de_expression.write("lrt <- glmLRT(fit, contrast=c{0})\n".format(tuple(contr)))
                    de_expression.write("FDR <- p.adjust(lrt$table$PValue, method='BH')\n")
                    de_expression.write("summary(dt <- decideTests(lrt))\n")
                    de_expression.write("isDE <- as.logical(dt)\n")
                    de_expression.write("DEnames <- rownames(y)[isDE]\n")
                    de_expression.write("outfile <- cbind(cpm(y), lrt$table, FDR)\n")
            return

    def runRscript(self, analysislocation):
        Rscriptloc = os.path.join(analysislocation, 'DEanalysis', 'DEexpression.r')
        subprocess.Popen("R --vanilla --slave < {0}".format(Rscriptloc), shell=True).wait()
        print "Analysis complete. Thank you for using SPARTA."
        return

    def de_analysis(self, analysislocation):
        self.removenoncountdata(analysislocation)
        conditions_list = self.getuserinput(analysislocation)
        self.writeRscript(analysislocation, conditions_list)
        self.runRscript(analysislocation)
        return

    def de_analysis_noninteractive(self, analysislocation, conditions_list):
        self.removenoncountdata(analysislocation)
        self.writeRscript(analysislocation, conditions_list)
        self.runRscript(analysislocation)
        return
