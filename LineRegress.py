#created by Brian Trethewey
#
#neGrapher is primary interface for graphing data
#neStats is primary interface for creating statistics output file for a dataset


import ConfigParser

from numpy import array, math
from scipy import stats
import matplotlib.pyplot as plt
from numpy import mean, median
import csv

#function to perform the linear regression and store the results in a dictionary
def lineRegress(linePoints):
    xArray, yArray = _pointsToVectors(linePoints)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xArray, yArray)
    result  = {}
    result["slope"] = slope
    result["intercept"] = intercept
    result["r_val"] = r_value
    result["p_val"] = p_value
    result["std_err"] = std_err
    return result

#helper function to create lines from slope, intercept data
def _getGraphLine(slope, intercept, xVctr):
    pointVctr = []
    for x in xVctr:
        pointVctr.append((x,slope*x + intercept))
    return pointVctr


def _pointsToVectors(points):
    xVctr = [point[0] for point in points]
    yVctr = [point[1] for point in points]
    return xVctr, yVctr


#function that implements the matplotlib and creates the graph object, also shows or saves the object to file
#dest parameter:  if show graph is printed to screen using .show() otherwise expects a filename with extension and saves the file in that format as per .savefig()
def createGraph(lineArray, title = None, xlab = None, yLab= None, colorVctr = None, styleVctr = None, dest = "show", xLim = None, yLim = None):
    lineCount = len(lineArray)
    if colorVctr and len(colorVctr) < lineCount:
        difference = lineCount-len(colorVctr)
        for  x in range(difference):
            colorVctr.append(colorVctr[len(colorVctr)-1])
    if styleVctr and len(styleVctr) < lineCount:
        difference = lineCount-len(styleVctr)
        for  x in range(difference):
            styleVctr.append(styleVctr[len(styleVctr)-1])

    if not colorVctr:
        colorFlag = ""
    if not styleVctr:
        styleFlag = ""
    for yVctrIdx in range(lineCount):
        if colorVctr:
            colorFlag = colorVctr[yVctrIdx]
        if styleVctr:
            styleFlag = styleVctr[yVctrIdx]
        argFlag = colorFlag+styleFlag
        xvctr, yvctr = _pointsToVectors(lineArray[yVctrIdx])

        xVctr = array(xVctr)
        yVctr = array(yvctr)
        plt.plot(xVctr,yVctr,argFlag)
    if title:
        plt.title(title)
    if xlab:
        plt.xlabel(xlab)
    if yLab:
        plt.ylabel(yLab)
    if xLim:
        plt.xlim(xLim)
    if yLim:
        plt.ylim(yLim)

    if dest == "show":
        plt.show()
    else:
        plt.savefig(dest, bbox_inches='tight')

#method to get teh confidence interval around the Slope of the regression
#uses the formula t((1-alpha/2):DoF)(s(b1))
#alpha: desired probability
#linePoints: list of touples defining the x and y coordinates of a point
#returns 3 variables, the slope of the regression, the intercept of the regression, and a touple containing the upper and lower bounds of the confidence interval of the slope
def slopeConfidence(alpha, linePoints):
    if len(linePoints)>2:
        return "Error: not enough points for calculation"
    if len(linePoints)==2:
        regression = lineRegress(linePoints)
        return regression["slope"], regression["intercept"],(regression["slope"],regression["slope"])
    #get linear regression for points
    regression = lineRegress(linePoints)
    #get Tscore
    tScore = stats.t.ppf(1-(alpha/2), len(linePoints)-2)


    #get s(b1)  == (MSE)/sum(xi-mean(x))^2)
    xVctr, yVctr = _pointsToVectors(linePoints)
    MSE = _MSE(regression["slope"], regression["intercept"],linePoints)
    xMean = mean(xVctr)
    xMeanDiffArray = []
    for x in xVctr:
        xDiff = x - xMean
        xDiffSq = xDiff *xDiff
        xMeanDiffArray.append(xDiffSq)
    xDiffSum = sum(xMeanDiffArray)
    sSqVal = MSE/xDiffSum
    sVal = math.sqrt(sSqVal)
    deltaConfidence = tScore*sVal
    confidenceInterval = (regression["slope"]-deltaConfidence,regression["slope"]+deltaConfidence)
    return regression["slope"], regression["intercept"], confidenceInterval


#get MSE of a regression
#slope: slope of calculated regression
#intercept: intercept of calculated regression
#yPoints: An array of the points used to create the regression
#xVctr: an array of the values
#returns MSE
def _MSE(slope, intercept,linePoints):
    errorArray = []

    #get  sigma error squared
    for point in range(len(linePoints)):
        xVal = point[0]
        yVal = point[1]
        expectedY = slope * xVal + intercept
        difference  = yVal - expectedY
        squareDifference  = difference*difference
        errorArray.append(squareDifference)
    errorSum = sum(errorArray)

    #devide by DoF (n-2)
    MSE = errorSum/(len(linePoints)-2)
    return MSE





#helper function to get means, extracted in case changes or medians are wanted instead
def _getExpectedLineValue(vctr):
    result = mean(vctr)
    return result

# xVctr:        List/Array of values representing the X axis(will be shared between all YPoints arrays)
# yPoints:      List/Array of Arrays composed of the Y values at each value in xVctr for each datastream
# expectedSlope:Optional Value for the slope of any Control or Expected line to compare to,
#       if None(default) no line produced
#       if value creates a line with that slope in red
#       if "auto" creates a line with slope =  average slope of all lines
def _NeRegressionGraphCalc(dataVctrs, expectedSlope = None):
    #get linear regression stats for all datasets
    LineStats = []
    for line in dataVctrs:
        data = lineRegress(line)
        LineStats.append(data)


    lineVctrs =[]
    colorVctr = []
    styleVctr = []

    #creates expected slope line for comparisons
    if expectedSlope:
        #get all slope and intercept values to get means
        slopes = []
        intercepts = []
        for statDict in LineStats:
            slopes.append(statDict["slope"])
            intercepts.append(statDict["intercept"])

        #get average intercept, this is important for placing the expected slope line
        expectedIntercept = _getExpectedLineValue(intercepts)

        #get auto generated expected SLope
        if expectedSlope == "auto":
            expectedSlope = _getExpectedLineValue(slopes)

        #make expected line for plotting
        #todo WTF do i do here??
        lineVctrs.append(x,_getGraphLine(expectedSlope, expectedIntercept, xVctr))
        colorVctr.append("r")
        styleVctr.append("-")

    for statDict in LineStats:
        slope = statDict["slope"]
        intercept = statDict["intercept"]
        #todo WTF do i do here??
        linePoints  = _getGraphLine(slope, intercept, xVctr)
        lineVctrs.append(linePoints)
        colorVctr.append("b")
        styleVctr.append("--")
    return lineVctrs, colorVctr,styleVctr

#combines linear regression and create graph into one function
def neGraphMaker(pointsVctrs, expectedSlope = None,title = None, xlab = None, yLab= None, dest = "show", xLim = None, yLim = None):
    lines, colors, styles = _NeRegressionGraphCalc(pointsVctrs, expectedSlope)
    createGraph(lines, colorVctr=colors, styleVctr=styles, title=title, xlab=xlab,yLab=yLab, dest=dest, xLim = xLim, yLim = yLim)

#reads in data fron neEst file outputs
def neFileRead(filename):
    fileBuffer = open(filename, "rb")
    replicateData = csv.DictReader(fileBuffer, delimiter="\t", quotechar="\"")
    dataDict = {}
    popNum = 0
    for item in replicateData:
        replicateNum = int(item['original_file'])
        popNum = int(item['pop'])
        neEst = float(item['est_ne'])
        if  not replicateNum in dataDict:
            dataDict[replicateNum] = {}
        dataDict[replicateNum][popNum] = neEst
    replicateKeys = dataDict.keys()
    resultTable = []
    for replicate in replicateKeys:
        replicateVctr = []
        replicateDict = dataDict[replicate]
        popKeys = replicateDict.keys()
        popKeys.sort()
        for popKey in popKeys:
            #print popKey
            replicateVctr.append((popKey,replicateDict[popKey]))
        resultTable.append(replicateVctr)
    return resultTable


#Method to read in a graph config file and return a dictionary of
def neConfigRead(filename):
    configDict = {}
    title =  None
    xLab = None
    yLab = None
    setExpected = None
    destination = "show"
    destType = "show"
    xLims =None
    yLims = None
    autoFlag = False
    config = ConfigParser.ConfigParser()
    config.readfp(open(filename))
    if config.has_section("labels"):
        if config.has_option("labels", "title"):
            title = config.get("labels", "title")
        if config.has_option("labels", "xLab"):
            xLab = config.get("labels", "xLab")
        if config.has_option("labels", "yLab"):
            yLab = config.get("labels", "yLab")
    if config.has_section("destination"):
        destType = config.get("destination","desttype")
        if destType != "show":
            destination = config.get("destination","destFile")
    if config.has_section("comparison"):
        if config.has_option("comparison", "auto"):
            autoFlag = config.getboolean("comparison", "auto")
            setExpected = "auto"
        if not autoFlag:
            if config.has_option("comparison", "lambda"):
                lambdaValue = config.getfloat("comparison", "lambda")
                setExpected = lambdaValue-1
            if config.has_option("comparison", "expectedSlope"):
                expectedSlope = config.getfloat("comparison", "expectedSlope")
                setExpected =  expectedSlope
    if config.has_section("limits"):
        if config.has_option("limits", "xMin") and config.has_option("limits", "xMax"):
            xMin = config.getfloat("limits", "xMin")
            xMax = config.getfloat("limits", "xMax")
            xLims = (xMin, xMax)
        if config.has_option("limits", "yMin")and config.has_option("limits", "yMax"):
            yMin = config.getfloat("limits", "yMin")
            yMax = config.getfloat("limits", "yMax")
            yLims =(yMin, yMax)
    configDict["title"]=title
    configDict["xLab"] = xLab
    configDict["yLab"] = yLab
    configDict["expected"] = setExpected
    configDict["dest"] = destination
    configDict["xLims"] = xLims
    configDict["yLims"] = yLims
    return configDict

#master function to create a graph from neEstimation data.
#neFile: filepath for the neEstimation output file desired.
#configFile filepath to configureation file containting parameteres for the graph (see example.cfg and example1.cfg,
#   this parameter and all feilds in the file are optional with what i considered the most relevant/base defaults)
def neGrapher(neFile, configFile):
    table = neFileRead(neFile)
    if not configFile:
        neGraphMaker(xVals,yTable)
        return True
    configs = neConfigRead(configFile)
    neGraphMaker(table,expectedSlope=configs["expected"],title= configs['title'],xlab=configs["xLab"],yLab=configs["yLab"],dest=configs["dest"],xLim=configs["xLims"],yLim=configs["yLims"])


#master function for creating a table of confidence intervals form neEstimation data
#neFile: filepath of neEstimation output file
#confidenceAlpha: level of confidence desired should be <0.5, most common are 0.05 and 0,1 for 95% and 90% respectivly
#outFileName: resulting file location for file results, will overwrite existing file.
#significantValue: value of comparison w/ regards to slope. should be 0 for every test, but can be changed if needed.
#testFlag: flag that disables file write and prints stats to console instead, used for test functions
def neStats(neFile,confidenceAlpha, outFileName = "neStatsOut.txt", significantValue = 0,testFlag = False):
    tableFormat = "{:<30}{:<30}{:<30}\n"
    tableString =tableFormat.format("Slope","Intercept","Confidence Interval")
    table = neFileRead(neFile)
    slopeVctr = []
    confidenceVctr = []

    for record in table:
        slope, intercept, confidence  = slopeConfidence(confidenceAlpha,record)
        tableString+=tableFormat.format(slope,intercept,confidence)
        slopeVctr.append(slope)
        confidenceVctr.append(confidence)
    maxSlope = max(slopeVctr)
    minSlope = min(slopeVctr)
    meanSlope = mean(slopeVctr)
    medSlope = median(slopeVctr)
    negativeCount=0
    zeroCount=0
    positiveCount=0
    for cI in confidenceVctr:
        if cI[0]>significantValue:
            positiveCount+=1
        elif cI[1]<significantValue:
            negativeCount+=1
        else:
            zeroCount+=1
    if testFlag:
        print slopeVctr
        print confidenceVctr

        print positiveCount
        print zeroCount
        print negativeCount
        return
    outFile = open(outFileName,"w")
    outFile.write("File:"+neFile+"\n")
    outFile.write("\n")
    outFile.write("Max Regression Slope:\t\t"+str(maxSlope)+"\n")
    outFile.write("Min Regression Slope:\t\t"+str(minSlope)+"\n")
    outFile.write("Mean Regression Slope:\t\t"+str(meanSlope)+"\n")
    outFile.write("Meadian Regression Slope:\t"+str(medSlope)+"\n")
    outFile.write("\n")
    outFile.write("Comparison to a slope of: "+str(significantValue))
    outFile.write("Positive Slopes:\t"+str(positiveCount)+"\t\tNeutral Slopes:\t"+str(zeroCount)+"\t\tNegative Slopes:\t"+str(negativeCount))
    outFile.write("\n\n")
    outFile.write(tableString)
    outFile.close()

if __name__ == "__main__":
    #Tests

    #Test Linear Regression
    print"lineRegress Tests"
    #Perfect Positive Regression
    xVct = range(10)
    yVct  = range(10)
    table = zip(xVct,yVct)
    result = lineRegress(table)
    assert result["slope"] == 1.0
    assert result["intercept"] == 0.0
    assert result["r_val"] ==1.0
    yVct = range(0,20,2)
    table = zip(xVct,yVct)
    result = lineRegress(xVct, yVct)
    assert result["slope"] == 2.0
    assert result["intercept"] == 0.0
    assert result["r_val"] == 1.0
    yVct = range(10,20,1)
    table = zip(xVct,yVct)
    result = lineRegress(xVct,yVct)
    assert result["slope"] == 1.0
    assert result["intercept"] == 10.0
    assert result["r_val"] ==1.0
    print "Perfect Positive regression passed"

    # Perfect Negative Regression
    xVct = range(10)
    yVct  = range(10,0,-1)
    table = zip(xVct,yVct)
    result = lineRegress(table)
    assert result["slope"] == -1.0
    assert result["intercept"] == 10.0
    assert result["r_val"] == -1.0
    yVct = range(20,0,-2)
    table = zip(xVct,yVct)
    result = lineRegress(table)
    assert result["slope"] == -2.0
    assert result["intercept"] == 20.0
    assert result["r_val"] ==-1.0
    yVct = range(0,-10,-1)
    table = zip(xVct,yVct)
    result = lineRegress(table)
    assert result["slope"] == -1.0
    assert result["intercept"] == 0.0
    assert result["r_val"] ==-1.0
    print "Perfect Negative regression passed"





    print "lineRegress Tests Passed"

    print "Confidence Tests "
    xVct = range(5)
    yVct  = range(5)
    table = zip(xVct,yVct)
    print slopeConfidence(0.05,table)
    xVct = range(5)
    yVct  = range(5)
    yVct[0]-=1
    table = zip(xVct,yVct)
    print slopeConfidence(0.01,table)
    print slopeConfidence(0.05, table)
    print slopeConfidence(0.1,table)
    xVct = range(10)
    yVct  = range(10)
    yVct[0]-=1
    table = zip(xVct,yVct)
    print slopeConfidence(0.01,table)
    print slopeConfidence(0.05, table)
    print slopeConfidence(0.1,table)
    print "Confidence Tests passed"

    print"getLineGraph Tests"
    xVct = range(10)

    assert _getGraphLine(0.0, 0.0, xVct) == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    print "slope 0 test passed"
    assert _getGraphLine(0.0, 6.0, xVct) == [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0]
    print "slope 0 intercept !=0 test passed"
    assert _getGraphLine(1.0, 0.0, xVct) == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    print"x=y case pass"
    assert _getGraphLine(1.0, 1.0, xVct) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    print "intercept != 0 test passed"
    assert _getGraphLine(4.0, 1.0, xVct) == [1.0, 5.0, 9.0, 13.0, 17.0, 21.0, 25.0, 29.0, 33.0, 37.0]
    print "slope 4 intercept 1 test passed"
    assert _getGraphLine(-1.0, 0.0, xVct) == [0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0]
    print "negative slope test passed"
    assert _getGraphLine(-2.0, 3.0, xVct) == [3.0, 1.0, -1.0 , -3.0, -5.0, -7.0, -9.0, -11.0, -13.0, -15.0]
    print "negative slope intercept !=0 test passed"
    assert _getGraphLine(2.0, -3.0, xVct) == [-3.0, -1.0, 1.0 , 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    print "negative intercept !=0 test passed"


    xVct = range(3,13,1)
    assert _getGraphLine(1.0, 0.0, xVct) == [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    print "non 0 positive x start passed"
    xVct = range(3,13,1)
    assert _getGraphLine(1.0, -3.0, xVct) == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    print "non 0 x start non 0 intercept passed"

    xVct = range(10,0,-1)
    assert _getGraphLine(1.0, 0.0, xVct) == [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
    print "negative step x axis passed"

    xVct = range(0, -10, -1)
    assert _getGraphLine(1.0, 0.0, xVct) == [0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0]
    print "negative step x axis passed"

    xVct = range(-7,3,1)
    assert _getGraphLine(1.0, 0.0, xVct) == [-7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
    print "negative start x axis passed"

    xVct = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    assert _getGraphLine(1.0, 0.0, xVct) == [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    print "non integer x axis passed"

    print "getGraphLine Passed"


    xVct = range(0,20,2)
    yVcts = []
    yVcts.append(_getGraphLine(-0.1, 0, xVct))
    yVcts.append(_getGraphLine(-0.005, 0, xVct))
    yVcts.append(_getGraphLine(-0.15, 10, xVct))
    yVcts.append(_getGraphLine(0.25, 10, xVct))


    print "createGraph Test"
    colorStyles = ["", ""]
    lineStyles = ["", ""]
    title = "testTitle"
    x = "TestX"
    y = "TestY"
    createGraph(xVct,yVcts,colorVctr=colorStyles, styleVctr=lineStyles, title=title, xlab=x, yLab=y)


    print "NeRegressionGraphMaker"

    testYs = [[0,0,0,0,0],[1,2,3,4,5],[-1,-2,-3,-4,-5],[1,3,1,3,1]]
    testX = [1,2,3,4,5]

    lines, color, styles = _NeRegressionGraphCalc(testX, testYs)

    createGraph(testX,lines,colorVctr=color,styleVctr=styles)

    testYs = [[0,0,0,0,0],[1,2,3,4,5],[-1,-2,-3,-4,-5],[1,3,0,2,-1]]
    lines, color, styles = _NeRegressionGraphCalc(testX, testYs, -0.333)

    createGraph(testX,lines,colorVctr=color,styleVctr=styles)

    testArray = [2,5,3,6,1,4]
    testArray.sort()
    print  testArray

    xVals, yTable = neFileRead("testData.txt")

    print xVals
    print yTable
    xVals = xVals[1:]
    temp = yTable[0]

    yTable = [temp[1:]]
    print xVals
    print yTable


    neGraphMaker(xVals,yTable)

    configwrite= ConfigParser.ConfigParser()
    configwrite.add_section("labels")
    configwrite.set("labels","title", "titletext")
    configwrite.set("labels", "xLab", "xLabel")
    configwrite.set("labels", "yLab", "yLabel")
    configwrite.add_section("destination")
    configwrite.set("destination","desttype", "PDF")
    configwrite.set("destination","destFile", "test.pdf")
    configwrite.add_section("comparison")
    configwrite.set("comparison", "auto", "True")
    configwrite.set("comparison", "expectedSlope", -0.1)
    configwrite.write(open("example1.cfg","w"))

    print neConfigRead("example1.cfg")
    print "test master methods"
    neGrapher("testData.txt","example1.cfg")
    neStats("testData.txt",0.1,testFlag=True)
    neStats("testData.txt",0.05)