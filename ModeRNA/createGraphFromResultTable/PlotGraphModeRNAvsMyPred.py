# script takes as input the output from evaluate_results.sh
# choose PlotGraph2 or PlotGraph1 method to be executed  


class Prediction(object):
    def __init__(self, target, template, length, similarity, gap, myRmsd, modeRnaRmsd):
        self.template = template
        self.target = target
        self.similarity = float(similarity)
        self.length = int(length)
        if modeRnaRmsd == "N/A":
            self.modeRnaRmsd = modeRnaRmsd
        else:
            self.modeRnaRmsd = float(modeRnaRmsd)
        if myRmsd == "N/A":
            self.myRmsd = myRmsd
        else:
            self.myRmsd = float(myRmsd)
        self.gap = float(gap)


class ListOfPredictions(object):
    predList = []
    def AppendPrediction(self, prediction):
        """
        :type prediction: Prediction
        """
        self.predList.append(prediction)
    def ComputeAverangeRMSDs(self, minLength, maxLength, maxRMSD, minSim, maxSim, minGap, maxGap):
        count = 0
        sumMy = 0
        countUnderRMSD = 0
        sumModeRNA = 0
        for p in self.predList:
            if (p.length < minLength) or (p.length > maxLength):
                continue
            if (p.similarity < minSim) or (p.similarity > maxSim):
                continue
            if (p.gap < minGap) or (p.gap > maxGap):
                continue
            if p.myRmsd < maxRMSD and p.myRmsd != "N/A":
                countUnderRMSD = countUnderRMSD + 1
            if (p.myRmsd == "N/A") or (p.myRmsd > maxRMSD) or (p.modeRnaRmsd == "N/A") or (p.modeRnaRmsd > maxRMSD):
                continue
            count = count + 1
            sumMy = sumMy + float(p.myRmsd)
            sumModeRNA = sumModeRNA + float(p.modeRnaRmsd)
        myAverage = sumMy / count
        modeRnaAverage = sumModeRNA / count
        return {"myAverage":myAverage, "modeRnaAverage":modeRnaAverage, "count":count, "curmsd": countUnderRMSD}


def LoadInputTable(tableName):
    inputTable = file(tableName, 'r')
    predictionList = ListOfPredictions()
    for line in inputTable: # expect lines like: TARGET  TEMPLATE    LEN    SIMIL      GAP     MY_RMSD MODERNA_RMSD
        lineSplited = line.split()
        if lineSplited[0] == "TARGET":
            continue
        predictionList.AppendPrediction(Prediction(lineSplited[0], lineSplited[1], lineSplited[2], lineSplited[3][:-1], lineSplited[4][:-1], lineSplited[5], lineSplited[6]))
    #ComputeAverangeRMSDs(self, minLength, maxLength, maxRMSD, minSim, maxSim, minGap, maxGap):
    result = predictionList.ComputeAverangeRMSDs(0, 9999, 30, 0, 100, 0, 100)
    print "myRMSD: " + str(result["myAverage"])
    print "modernaRMSD: " + str(result["modeRnaAverage"])
    print "count: " + str(result["count"])
    print "count under certain RMSD: " + str(result["curmsd"])

    result = predictionList.ComputeAverangeRMSDs(50, 100, 30, 0, 100, 0, 100)
    print "myRMSD: " + str(result["myAverage"])
    print "modernaRMSD: " + str(result["modeRnaAverage"])
    print "count: " + str(result["count"])
    print "count under certain RMSD: " + str(result["curmsd"])
    result = predictionList.ComputeAverangeRMSDs(101, 500, 30, 0, 100, 0, 100)
    print "myRMSD: " + str(result["myAverage"])
    print "modernaRMSD: " + str(result["modeRnaAverage"])
    print "count: " + str(result["count"])
    print "count under certain RMSD: " + str(result["curmsd"])

    return predictionList

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter



def PlotGraph2(predictionList):
    # data
    x = []
    y = []
    color = []
    for a in predictionList.predList:
        if (a.myRmsd != "N/A") and (a.myRmsd < 30):
            y.append(a.myRmsd)
            x.append(a.length)
            if a.similarity <= 75 :
                color.append('green')  # sim is < 75 => green point
            else:
                color.append('red')  # sim is > 75 =>   red point

    nullfmt = NullFormatter()  # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    #rect_histx = [left, bottom_h, width, 0.2]
    #rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)

    plt.xlabel("Length")
    plt.ylabel("RMSD")

    #axHistx = plt.axes(rect_histx)
    #axHisty = plt.axes(rect_histy)

    # no labels
    #axHistx.xaxis.set_major_formatter(nullfmt)
    #axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, color=color)

    # now determine nice limits by hand:
    #binwidth = 0.25
    #xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
    #lim = (int(xymax / binwidth) + 1) * binwidth

    #axScatter.set_xlim((0, lim))
    axScatter.set_ylim((0, 30))

    #bins = np.arange(0, lim + binwidth, binwidth)
    #axHistx.hist(x, bins=bins)
    #axHisty.hist(y, bins=bins, orientation='horizontal')

    #axHistx.set_xlim(axScatter.get_xlim())
    #axHisty.set_ylim(axScatter.get_ylim())

    plt.show()


def PlotGraph1(predictionList):
    # data
    x = []
    y = []
    color = []
    for a in predictionList.predList:
        if (a.myRmsd != "N/A") and (a.modeRnaRmsd != "N/A") and (a.myRmsd < 30 ) and (a.modeRnaRmsd < 30):
            x.append(a.myRmsd)
            y.append(a.modeRnaRmsd)
            if a.length < 100:
                color.append('green') # len is < 50 => green point
            else:
                color.append('red')   # len is > 50 =>   red point

    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02


    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)

    plt.xlabel("Ours")
    plt.ylabel("ModeRNA")

    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, color=color)

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
    lim = (int(xymax/binwidth) + 1) * binwidth

    axScatter.set_xlim((0, lim))
    axScatter.set_ylim((0, lim))

    bins = np.arange(0, lim + binwidth, binwidth)
    axHistx.hist(x, bins=bins)
    axHisty.hist(y, bins=bins, orientation='horizontal')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    plt.show()

def CompareResults(list_w_o_sec_str, list_with_sec_str):
    # idea: put all records of second list into dictionary and them index them by the first list - and compare
    dictionary_with_sec_str = {}
    for record in list_with_sec_str.predList:
        if record.myRmsd != 'N/A' and record.myRmsd < 50:
            dictionary_with_sec_str[record.target + record.template] = record.myRmsd
    for record in list_w_o_sec_str.predList:
        if record.myRmsd != 'N/A' and record.myRmsd < 50:
            if dictionary_with_sec_str.has_key(record.target + record.template):
                if dictionary_with_sec_str[record.target + record.template] > record.myRmsd + 3:
                    print  'Target: ' + record.target + ', Tempalte: ' + record.template + ', RMSD_w_o_sec_str: ' + str(record.myRmsd) + ', RMSD_with_sec_str: ' + str(dictionary_with_sec_str[record.target + record.template])


predictionList = LoadInputTable("result_table.txt")
predictionList_w_o_sec_str = LoadInputTable("result_table_w_o_sec_str.txt")
CompareResults(predictionList_w_o_sec_str, predictionList)
#PlotGraph1(predictionList)
#PlotGraph2(predictionList)