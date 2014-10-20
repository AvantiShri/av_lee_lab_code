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
    , SIGNED_LOG_PVAL_METASCORES = 'signedLogPvalMetascores'
    , GENE_ID_COLUMN = 'geneIdColumn'
    );
VALUE_COLUMNS_TO_STORE_ATTR = util.enum(
    VALUE_NAME = 'valueName'
    , COLUMN_NAME = 'columnName'
    , TYPE = 'type'
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
    MATURATION_SCORE = "maturationScore"
    ,D7_V_SHAM_SCORE = "d7RvsShamScore"
    ,DEDIFF_SCORE = "dediffScore"
    ,ESC_DIFF_SCORE = "escDiffScore" 
);


def main():
    parser = argparse.ArgumentParser(description="read in the lee lab files in preparation for statistical tests");
    parser.add_argument("--inputConfigs", nargs="+", required=True);
    args = parser.parse_args();
    
    genes = processInput(args).values();

    numGenesToAverage = 500;
    numIterations = 10000;

    compareToScore(genes, SCORE_NAMES.ESC_DIFF_SCORE
        , [Settings(scoreToCompare=SCORE_NAMES.D7_V_SHAM_SCORE, reverseScoreToCompareOrder=True, greaterThanThreshold=False)]
        , numGenesToAverage, numIterations);
    
    compareToScore(genes, SCORE_NAMES.MATURATION_SCORE
        , [Settings(scoreToCompare=SCORE_NAMES.DEDIFF_SCORE, reverseScoreToCompareOrder=True, greaterThanThreshold=False)]
        , numGenesToAverage, numIterations);


class Settings(object):
    def __init__(self, scoreToCompare, reverseScoreToCompareOrder, greaterThanThreshold):
        self.scoreToCompare = scoreToCompare;
        self.reverseScoreToCompareOrder = reverseScoreToCompareOrder;
        self.greaterThanThreshold = greaterThanThreshold;

def compareToScore(genes, scoreToCompareTo, settings, numGenesToAverage, numIterations):
    (genes, monteCarloDist) = filterAndMonteCarlo(genes, scoreToCompareTo, numGenesToAverage, numIterations);
    for setting in settings:
        score = twoScoreComparison(
            genes = genes
            , scoreToCompare=setting.scoreToCompare
            , reverseScoreToCompareOrder=setting.reverseScoreToCompareOrder
            , scoreToCompareTo=scoreToCompareTo
            , numToAverage = numGenesToAverage
            , monteCarloDistribution = monteCarloDist
            , greaterThanThreshold = setting.greaterThanThreshold
            );    
    

def filterAndMonteCarlo(genes, scoreName, numGenesToAverage, numIterations):
    genes = filterOutZeros(genes, scoreName);
    monteCarloDist = getMonteCarloDist(genes, scoreName, numGenesToAverage, numIterations);
    return (genes, monteCarloDist);
    

def getMonteCarloDist(genes,scoreName,numGenesToSample,numIterations):
    monteCarloDist = monteCarloDistribution(genes, scoreName, numGenesToSample, numIterations);
    return monteCarloDist;

def twoScoreComparison(
            genes
            , scoreToCompare
            , scoreToCompareTo
            , numToAverage
            , monteCarloDistribution
            , reverseScoreToCompareOrder=False
            , greaterThanThreshold=True
    ):
    genes = filterOutZeros(genes, scoreToCompareTo);
    genes = sortByScore(genes, scoreToCompare, reverseScoreToCompareOrder);
    threshold = averageTopN(genes, scoreToCompareTo, numToAverage);
    zScore = getZscore(monteCarloDistribution,threshold);
    proportionGreaterThan = getProportionComparedTo(monteCarloDistribution, threshold, greaterThanThreshold);
    
    fileName = scoreToCompare+"_vs_"+scoreToCompareTo+"_"+("descen" if reverseScoreToCompareOrder else "ascen");
 
    score = stats.TestResult(proportionGreaterThan, "Monte carlo with "+str(len(monteCarloDistribution))+" samples",testStatistic=zScore);  
    print fileName+" "+str(score);    
    plt.hist(monteCarloDistribution, bins=np.linspace(min(monteCarloDistribution),max(monteCarloDistribution),100));
    plt.plot([threshold,threshold],[0,100]);
    #plt.show();
    plt.savefig(fileName+".svg");
    plt.savefig(fileName+".png");
    plt.close("all");

def filterOutZeros(genes, scoreName):
    return [x for x in genes if x.getAttribute(scoreName) != 0];

def sortByScore(genes, scoreName, reverse=False):
    genes.sort(reverse=reverse,key=lambda x: x.getAttribute(scoreName));
    return genes;

def averageTopN(genes, scoreName, numToAverage):
    theSum = 0;
    for i in range(numToAverage):
        theSum += genes[i].getAttribute(scoreName);
    return float(theSum)/numToAverage;

def monteCarloDistribution(genes,scoreName, numGenesToSample, numIterations):
    def monteCarloAction():
        sample = random.sample(genes, numGenesToSample);
        return averageTopN(sample, scoreName, numGenesToSample);
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
        geneIdColumn = configJson[CONFIG_ATTR.GENE_ID_COLUMN];    
        valueColumnsToStore = configJson[CONFIG_ATTR.VALUE_COLUMNS_TO_STORE];
        signedLogPvalMetascores = configJson[CONFIG_ATTR.SIGNED_LOG_PVAL_METASCORES];

        def actionOnDictionary(inp):
            geneId = inp[geneIdColumn];
            if geneId not in geneDictionary:
                geneDictionary[geneId] = util.Entity(geneId);
            theGene = geneDictionary[geneId];
            for valueColumnToStore in valueColumnsToStore:
                theType = valueColumnToStore[VALUE_COLUMNS_TO_STORE_ATTR.TYPE] if VALUE_COLUMNS_TO_STORE_ATTR.TYPE in valueColumnToStore else None;
                columnName = valueColumnToStore[VALUE_COLUMNS_TO_STORE_ATTR.COLUMN_NAME];
                theVal = inp[columnName];
                valueName = valueColumnToStore[VALUE_COLUMNS_TO_STORE_ATTR.VALUE_NAME] if VALUE_COLUMNS_TO_STORE_ATTR.VALUE_NAME in valueColumnToStore else columnName;
                theGene.addAttribute(
                    valueName
                    , theVal if theType is None else util.transformType(theVal, theType)
                );
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

                    signedLogPval = theGene.getAttribute(pvalName);
                    signedLogPval = abs(math.log(signedLogPval) if fcIsLog == False else signedLogPval);
                    fcThreshold = 0 if fcIsLog else -1;
                    signedLogPval = signedLogPval if theGene.getAttribute(fcColName) >= fcThreshold else -1*signedLogPval;
                    signedLogPval = signedLogPval if invertSign == False else -1*signedLogPval;
                    metascore += signedLogPval;
                theGene.addAttribute(metascoreName, metascore);
        
        def actionFromTitle(title):
            dictionaryFromLine = fp.lambdaMaker_dictionaryFromLine(title);
            def theAction(x,i):
                actionOnDictionary(dictionaryFromLine(x));
            return theAction;

        fp.performActionOnEachLineOfFile(fileHandle=fp.getFileHandle(inputFile), actionFromTitle=actionFromTitle, ignoreInputTitle=True);

    return geneDictionary;

main();
