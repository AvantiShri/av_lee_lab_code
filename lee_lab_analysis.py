#!/usr/bin/python
import sys;
import os;
scriptsDir = os.environ.get("UTIL_SCRIPTS_DIR");
if (scriptsDir is None):
	raise Exception("Please set environment variable UTIL_SCRIPTS_DIR");
sys.path.insert(0,scriptsDir);
import pathSetter;
import argparse;
import util;
import fileProcessing as fp;
import random;
import stats;
import matplotlib.pyplot as plt;
import numpy as np;
import math;
from scipy.stats import shapiro;

#Next:
#Extract the values from the gene dictionary

#set up an enum for the names of the esc attribute and stuff.

#processing pipeline:
#filter out values where score is zero for esc score
#sort by d7 score
#

CONFIG_ATTR = util.enum(
    INPUT_FILE = 'inputFile'
    , VALUE_COLUMNS_TO_STORE = 'valueColumnsToStore'
    , GENE_ID_COLUMN = 'geneIdColumn'
    , FILE_TYPE = 'fileType'
    , SAMPLE1COL = 'sample_1'
    , SAMPLE2COL = 'sample_2'
    , RENAMINGS = 'renamings' 
    );
METASCORES_ATTR = util.enum(
    SIGNED_LOG_PVAL_METASCORES = 'signedLogPvalMetascores'
);
VALUE_COLUMNS_TO_STORE_ATTR = util.enum(
    VALUE_NAME = 'valueName'
    , COLUMN_NAME = 'columnName'
    , TYPE = 'type'
    , FLIP_SIGN_ON_INVERT = 'flipSignOnInvert' 
);
SIGNED_LOG_PVAL_METASCORES_ATTR = util.enum(
    SCORE_NAME = 'scoreName'
    , PVAL_FOLDCHANGE_PAIRS = 'pvalFoldChangePairs'
);
PVAL_FOLD_CHANGE_PAIRS_ATTR = util.enum(
    PVAL_NAME = 'pvalName'
    , FC_COL_NAME = 'fcColName'
    , INVERT_SIGN = 'invertSign'
    , PVAL_IS_LOG = 'pvalIsLog' #default to false
    , FC_IS_LOG = 'fcIsLog' #default to true
);
SCORE_NAMES = util.enum(
    MATURATION_SCORE_OLDADULT = "maturationScoreOldAdult"
    ,MATURATION_SCORE_NEWADULT = "maturationScoreNewAdult"
    ,P0P4P7_SCORE = "P0P4P7Score"
    ,P0_V_P7_SCORE = "P0vsP7Score"
    ,P0_V_ADULT_SCORE = "P0vsAdultScore"
    ,P7_V_ADULT_SCORE = "P7vsAdultScore"
    ,D7_V_SHAM_SCORE = "d7RvsShamScore"
    ,DEDIFF_SCORE = "dediffScore"
    ,ESC_DIFF_SCORE = "escDiffScore"
    ,ESC_V_CP = "ESCvsCP"
    ,ESC_V_CM = "ESCvsCM" 
);

zeroPvalSubstitute = 10**-250;

def main():
    parser = argparse.ArgumentParser(description="read in the lee lab files in preparation for statistical tests");
    parser.add_argument("--inputConfigs", nargs="+", required=True);
    parser.add_argument("--metascoresConfig", required=True);
    parser.add_argument("--prefix",required=True);
    args = parser.parse_args();
    
    genes = processInput(args).values();

    numGenesToAverage = 500;
    numIterations = 10000;

    d7vShamSettings = Settings(SCORE_NAMES.D7_V_SHAM_SCORE, args.prefix, 
                [Setting(reverseScoreToCompareOrder=True, greaterThanThreshold=False)
                , Setting(reverseScoreToCompareOrder=False, greaterThanThreshold=True)]);

    outputFile = args.prefix+"_scores.tsv";
    scoresToPrint = [SCORE_NAMES.D7_V_SHAM_SCORE, SCORE_NAMES.ESC_DIFF_SCORE, SCORE_NAMES.MATURATION_SCORE_OLDADULT, SCORE_NAMES.MATURATION_SCORE_NEWADULT, SCORE_NAMES.DEDIFF_SCORE]
    #for score in scoresToPrint:
    #    computePercentiles(genes,score);
    #util.printAttributes(genes, scoresToPrint+[x+"_perc" for x in scoresToPrint], outputFile);
    
    scoresToTestWith = [SCORE_NAMES.ESC_DIFF_SCORE, SCORE_NAMES.MATURATION_SCORE_OLDADULT, SCORE_NAMES.MATURATION_SCORE_NEWADULT, SCORE_NAMES.DEDIFF_SCORE];
    #scoresToTestWith =  scoresToTestWith+[x+"_perc" for x in scoresToTestWith];
    for scoreToTestWith in scoresToTestWith:
        compareToScore(genes, scoreToTestWith, d7vShamSettings, numGenesToAverage, numIterations);


def computePercentiles(genes, scoreName):
    genes = sortByScore(genes, scoreName);
    i = 0;
    while i < len(genes):
        genes[i].addAttribute(scoreName+"_perc", (100*float(i+1))/len(genes));
        i += 1;

class Settings(object):
    def __init__(self, scoreToCompare, prefix, settingsArr):
        self.scoreToCompare = scoreToCompare;
        self.prefix = prefix;
        self.settingsArr = settingsArr;

class Setting(object):
    def __init__(self, reverseScoreToCompareOrder, greaterThanThreshold):
        self.reverseScoreToCompareOrder = reverseScoreToCompareOrder;
        self.greaterThanThreshold = greaterThanThreshold;

def compareToScore(genes, scoreToCompareTo, settings, numGenesToAverage, numIterations):
    (genes, monteCarloDist) = filterAndMonteCarlo(genes,scoreToCompareTo,scoreToCompareTo,numGenesToAverage,numIterations);
    #(genes,monteCarloDist) = filterAndMonteCarlo(genes, scoreToCompareTo, settings.scoreToCompare, numGenesToAverage, numIterations);
    #print shapiro(monteCarloDist);
    for setting in settings.settingsArr:
        score = twoScoreComparison(
            genes
            , settings.prefix
            , scoreToCompare=settings.scoreToCompare
            , reverseScoreToCompareOrder=setting.reverseScoreToCompareOrder
            , scoreToCompareTo=scoreToCompareTo
            , numToAverage = numGenesToAverage
            , monteCarloDistribution = monteCarloDist
            , greaterThanThreshold = setting.greaterThanThreshold
            );    
    

def filterAndMonteCarlo(genes, scoreName, scoreToFilterZeros, numGenesToAverage, numIterations):
    genes = filterOutZeros(genes, scoreToFilterZeros);
    monteCarloDist = getMonteCarloDist(genes, scoreName, numGenesToAverage, numIterations);
    return (genes, monteCarloDist);
    

def getMonteCarloDist(genes,scoreName,numGenesToSample,numIterations):
    monteCarloDist = monteCarloDistribution(genes, scoreName, numGenesToSample, numIterations);
    return monteCarloDist;

def twoScoreComparison(
            genes
            , prefix
            , scoreToCompare
            , scoreToCompareTo
            , numToAverage
            , monteCarloDistribution
            , reverseScoreToCompareOrder=False
            , greaterThanThreshold=True
    ):
    #genes = filterOutZeros(genes, scoreToCompareTo);
    genes = sortByScore(genes, scoreToCompare, reverseScoreToCompareOrder);
    threshold = medianTopN(genes, scoreToCompareTo, numToAverage);
    zScore = getZscore(monteCarloDistribution,threshold);
    proportionGreaterThan = getProportionComparedTo(monteCarloDistribution, threshold, greaterThanThreshold);
    
    fileName = prefix+"_"+scoreToCompare+"_vs_"+scoreToCompareTo+"_"+("descen" if reverseScoreToCompareOrder else "ascen");
 
    score = stats.TestResult(proportionGreaterThan, "Monte carlo with "+str(len(monteCarloDistribution))+" samples",testStatistic=zScore, testContext="threshold: "+str(threshold));  
    plt.hist(monteCarloDistribution, bins=np.linspace(min(monteCarloDistribution),max(monteCarloDistribution),100));
    plt.plot([threshold,threshold],[0,1000]);
    #plt.show();
    plt.savefig(fileName+".svg");
    plt.savefig(fileName+".png");
    plt.close("all");
    print fileName+": "+str(score.testStatistic);
    print fileName+": "+str(score);

def filterOutZeros(genes, scoreName):
    return [x for x in genes if x.getAttribute(scoreName) != 0];

def sortByScore(genes, scoreName, reverse=False):
    genes.sort(reverse=reverse,key=lambda x: x.getAttribute(scoreName));
    return genes;

def medianTopN(genes, scoreName, numToGetMedOf):
    return np.median([x.getAttribute(scoreName) for x in genes[0:numToGetMedOf]]);

def monteCarloDistribution(genes,scoreName, numGenesToSample, numIterations):
    def monteCarloAction():
        sample = random.sample(genes, numGenesToSample);
        return medianTopN(sample, scoreName, numGenesToSample);
    return stats.monteCarlo(monteCarloAction, numIterations);

def getZscore(arr, val):
    avg = stats.mean(arr);
    sdev = stats.sdev(arr);
    return float(val - avg)/sdev;

def getProportionComparedTo(arr, val, greaterThan=True):
    numGE = 0;
    for elem in arr:
        if ((greaterThan and elem >= val) or ((greaterThan == False) and elem <= val)):
            numGE += 1;
    return float(numGE)/len(arr);


#write functions to:
# filter out all values where score was 0.
# sort by a given score
# take average of top n of a particular score in an array
# randomly sample n scores from array m times and generate the values of the distribution.
# compare a threshold to the values from a distribution to form a 'z' score, maybe also plot.
# copy top n of a particular score into an array (for wilcoxon comparison).
# wilcoxon comparison.

#...and that should cover it, kindof...

def processInput(args):
    geneDictionary = {};
    for inputConfig in args.inputConfigs:
        configJson = util.parseJsonFile(inputConfig);
        inputFile = configJson[CONFIG_ATTR.INPUT_FILE];
        fileType = configJson[CONFIG_ATTR.FILE_TYPE] if CONFIG_ATTR.FILE_TYPE in configJson else "default";
        geneIdColumn = configJson[CONFIG_ATTR.GENE_ID_COLUMN];    
        valueColumnsToStore = configJson[CONFIG_ATTR.VALUE_COLUMNS_TO_STORE];
        def actionOnDictionary(inp):
            geneId = inp[geneIdColumn];
            if geneId not in geneDictionary:
                geneDictionary[geneId] = util.Entity(geneId);
            theGene = geneDictionary[geneId];
            if (fileType == "default"):
                for valueColumnToStore in valueColumnsToStore:
                    (theVal, valueName) = basicValueColumnExtraction(valueColumnToStore, inp); 
                    theGene.addAttribute(valueName, theVal);
            elif (fileType == "cuffdiff"):
                sample1 = inp[CONFIG_ATTR.SAMPLE1COL];
                sample2 = inp[CONFIG_ATTR.SAMPLE2COL];
                renamings = configJson[CONFIG_ATTR.RENAMINGS] if CONFIG_ATTR.RENAMINGS in configJson else {};
                sample1 = renamings[sample1] if sample1 in renamings else sample1;
                sample2 = renamings[sample2] if sample2 in renamings else sample2;
                for valueColumnToStore in valueColumnsToStore:
                    (theVal, valueNameCore) = basicValueColumnExtraction(valueColumnToStore, inp);
                    valueName = valueNameCore+"_"+sample1+"_"+sample2;
                    theGene.addAttribute(
                        valueName, theVal
                    );
                    flipSignOnInvert = valueColumnToStore[VALUE_COLUMNS_TO_STORE_ATTR.FLIP_SIGN_ON_INVERT] if VALUE_COLUMNS_TO_STORE_ATTR.FLIP_SIGN_ON_INVERT in valueColumnToStore else False;
                    theVal = -1*theVal if flipSignOnInvert else theVal;
                    valueName = valueNameCore+"_"+sample2+"_"+sample1;
                    theGene.addAttribute(
                        valueName, theVal
                    );                        
            else:
                raise ValueError("Unrecognized fileType "+fileType); 
        
        def actionFromTitle(title):
            dictionaryFromLine = fp.lambdaMaker_dictionaryFromLine(title);
            def theAction(x,i):
                actionOnDictionary(dictionaryFromLine(x));
            return theAction;

        fp.performActionOnEachLineOfFile(fileHandle=fp.getFileHandle(inputFile), actionFromTitle=actionFromTitle, ignoreInputTitle=True);

    metascoresJson = util.parseJsonFile(args.metascoresConfig);
    signedLogPvalMetascores = metascoresJson[METASCORES_ATTR.SIGNED_LOG_PVAL_METASCORES];

    for theGene in geneDictionary.values():
        for signedLogPvalMetascore in signedLogPvalMetascores:
            metascoreName = signedLogPvalMetascore[SIGNED_LOG_PVAL_METASCORES_ATTR.SCORE_NAME];
            pvalFoldChangePairs = signedLogPvalMetascore[SIGNED_LOG_PVAL_METASCORES_ATTR.PVAL_FOLDCHANGE_PAIRS];
            metascore = 0;
            for pvalFoldChangePair in pvalFoldChangePairs:
                pvalName = pvalFoldChangePair[PVAL_FOLD_CHANGE_PAIRS_ATTR.PVAL_NAME];
                fcColName = pvalFoldChangePair[PVAL_FOLD_CHANGE_PAIRS_ATTR.FC_COL_NAME];
                invertSign = pvalFoldChangePair[PVAL_FOLD_CHANGE_PAIRS_ATTR.INVERT_SIGN] if PVAL_FOLD_CHANGE_PAIRS_ATTR.INVERT_SIGN in pvalFoldChangePair else False;
                pvalIsLog = pvalFoldChangePair[PVAL_FOLD_CHANGE_PAIRS_ATTR.PVAL_IS_LOG] if PVAL_FOLD_CHANGE_PAIRS_ATTR.PVAL_IS_LOG in pvalFoldChangePair else False;
                fcIsLog = pvalFoldChangePair[PVAL_FOLD_CHANGE_PAIRS_ATTR.FC_IS_LOG] if PVAL_FOLD_CHANGE_PAIRS_ATTR.FC_IS_LOG in pvalFoldChangePair else True; 

                if (theGene.hasAttribute(pvalName) == False):
                    next; #skip the iteration; pval insignificant, or something.
                else:
                    signedLogPval = theGene.getAttribute(pvalName);
                    signedLogPval = abs(math.log(zeroPvalSubstitute if signedLogPval < zeroPvalSubstitute else signedLogPval) if pvalIsLog == False else signedLogPval);
                    fcThreshold = 0 if fcIsLog else 1;
                    signedLogPval = signedLogPval if theGene.getAttribute(fcColName) >= fcThreshold else -1*signedLogPval;
                    signedLogPval = signedLogPval if invertSign == False else -1*signedLogPval;
                    metascore += signedLogPval;
            theGene.addAttribute(metascoreName, metascore);

    return geneDictionary;

def basicValueColumnExtraction(valueColumnToStore, inp):
    theType = valueColumnToStore[VALUE_COLUMNS_TO_STORE_ATTR.TYPE] if VALUE_COLUMNS_TO_STORE_ATTR.TYPE in valueColumnToStore else None;
    columnName = valueColumnToStore[VALUE_COLUMNS_TO_STORE_ATTR.COLUMN_NAME];
    theVal = inp[columnName];
    valueName = valueColumnToStore[VALUE_COLUMNS_TO_STORE_ATTR.VALUE_NAME] if VALUE_COLUMNS_TO_STORE_ATTR.VALUE_NAME in valueColumnToStore else columnName;    
    theVal = theVal if theType is None else util.transformType(theVal, theType)

    return (theVal, valueName); 

main();
