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
    def __init__(self):
        self.predList = []
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
    print "********ALL*******"
    print "myRMSD: " + str(result["myAverage"])
    print "modernaRMSD: " + str(result["modeRnaAverage"])
    print "count: " + str(result["count"])
    print "count under certain RMSD: " + str(result["curmsd"])

    result = predictionList.ComputeAverangeRMSDs(50, 100, 30, 0, 100, 0, 100)
    print "********50-100*******"
    print "myRMSD: " + str(result["myAverage"])
    print "modernaRMSD: " + str(result["modeRnaAverage"])
    print "count: " + str(result["count"])
    print "count under certain RMSD: " + str(result["curmsd"])
    result = predictionList.ComputeAverangeRMSDs(101, 500, 30, 0, 100, 0, 100)
    print "********101-500*******"
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
    # dictionary: {with secondary structure, without secondary structure, ModeRNA prediction}
    wss_woss_mrna_dict = {}
    for record in list_w_o_sec_str.predList:
        if record.myRmsd != 'N/A' and record.myRmsd < 50:
            if dictionary_with_sec_str.has_key(record.target + record.template):
                if record.modeRnaRmsd != 'N/A' and record.length < 1000: # ToDo: parametrizuj podla dlzky
                    wss_woss_mrna_dict[record.target + record.template] = [dictionary_with_sec_str[record.target + record.template] ,record.myRmsd, record.modeRnaRmsd, record.length]
                if dictionary_with_sec_str[record.target + record.template] > record.myRmsd + 3:
                    pass
                    #print 'Target: ' + record.target + ', Tempalte: ' + record.template + ', RMSD_w_o_sec_str: ' + str(record.myRmsd) + ', RMSD_with_sec_str: ' + str(dictionary_with_sec_str[record.target + record.template])
    return wss_woss_mrna_dict

def CompareAvarageResults(dict):
    print sum(dict[l][0] for l in dict) / len(dict)
    print sum(dict[l][1] for l in dict) / len(dict)
    print sum(dict[l][2] for l in dict) / len(dict)

# returns list: [RMSD_WITH_SECSTRU, RMSD_WO_SECSTRU, # of FalsePositives in secstru prediction divided by structure length]
def EnrichDictionaryWithSecStruFalsePositives(dictToEnrich, enrichByDict):
    enrichedList = []
    for key in dictToEnrich:
        # get target structure name
        targetName = key[:6]
        # select false positives and add it to corresponding line in dictionary
        try:
            enrichedList.append([dictToEnrich[key][0]
                                    , dictToEnrich[key][1]
                                    , float(enrichByDict[targetName][2])
                                    , dictToEnrich[key][3]])
            #/dictToEnrich[key][3]]) # put the division to third column and let the #ofFalsePositives be parametrized by length
        except:
            pass
            # it can happen: I do not have to have secStru of all records, template secstru is enough
    return enrichedList

# as input takes list of [RMSD_WITH_SECSTRU, RMSD_WO_SECSTRU, # of FalsePositives in secstru prediction]
def PlotGraptCompareWithAndWoSecStr(list):
    font = {'family': 'normal',
            'size': 25}
    plt.rc('font', **font)
    # count the diference between with and wo secondary struct prediction
    # negative number means that we went worse, positive that we went better
    # x axis will be # of false positived in secstru prediction
    yaxis = []
    xaxis = []
    label = []
    for rec in list:
        yaxis.append(float(rec[1] - rec[0]))
        xaxis.append(float(rec[2]))
        label.append(int(rec[3]))
    color = ['red' if l < 100 else 'green' for l in label]
    # linear regression
    x = np.array(xaxis)
    y = np.array(yaxis)
    fit = np.polyfit(x, y, 1)
    fit_fn = np.poly1d(fit)
    plt.plot(x, y, 'ro', x, fit_fn(x), '--k') #'ro',
    #plt.scatter(x, y, color=color)
    #plt.plot(x, fit_fn(x), '--k')
    plt.axhline(y = 0)
    plt.xlabel("FalsePositives in secondary structure prediction (# of basepairs)")
    plt.ylabel("Differnce between predicted structures (RMSD)")
    plt.show()


# returns dictionary like: structureName : TruePos FalseNeg FalsePos
# as input takes output file from CompareSecondaryStructs.py
def LoadSecstruComparationResults(filename):
    file = open(filename, 'rU')
    firstLine = True
    # structure TruePos FalseNeg FalsePos TrueNeg
    dict = {}
    for line in file:
        if firstLine:
            firstLine = False
            continue
        splitedLine = str.split(line)
        dict[splitedLine[0]] = [splitedLine[1], splitedLine[2], splitedLine[3]]
    file.close()
    return dict

def ComparationRelatedToSecStru():
    print '###################WITH SECONDARY STRUCTURE###########################'
    predictionList = LoadInputTable("result_table.txt")
    print '##################WITHOUT SECONDARY STRUCTURE############################'
    predictionList_w_o_sec_str = LoadInputTable("result_table_w_o_sec_str.txt")
    print '##############################################'
    dict = CompareResults(predictionList_w_o_sec_str, predictionList)
    print '##############################################'
    #CompareAvarageResults(dict)
    secStruComparation = LoadSecstruComparationResults("SecStruComparation.txt")
    enrichedList = EnrichDictionaryWithSecStruFalsePositives(dict, secStruComparation)
    PlotGraptCompareWithAndWoSecStr(enrichedList)

ComparationRelatedToSecStru()

#PlotGraph1(predictionList)
#PlotGraph2(predictionList)