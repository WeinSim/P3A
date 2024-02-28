# counts the number of points that lie between (slope * x + yInter - unc)
# and (slope * x + yInter + unc)
def numPointsInside(xvals, yvals, slope, yInter, unc):
    ret = 0
    for i in range(len(xvals)):
        if abs(xvals[i] * slope + yInter - yvals[i]) < unc:
            ret += 1
    return ret

# Calculates the uncertainty of the slope and y intercept by finding
# an interval which includes two thirds of the points. This is the uncertainty
# of the y intercept. The uncertainty of the slope is obtained
# by calculating the slopes of the diagonals of this uncertainty strip.
# Returns a tuple of the form (slope uncertainty, y inter uncertainty).
def linRegUncertainty(xvals, yvals, coefs):
    slope = coefs[0]
    yInter = coefs[1]
    numSteps = 10
    targetPercent = 68
    targetPoints = int(len(xvals) * targetPercent * 0.01)
    maxDist = 0 
    for i in range(len(xvals)):
        x = xvals[i]
        y = yvals[i]
        maxDist = max(maxDist, abs(x * slope + yInter - y))
    minUnc = 0
    maxUnc = maxDist + 0.1
    for i in range(numSteps):
        midUnc = (minUnc + maxUnc) / 2
        pointsInside = numPointsInside(xvals, yvals, slope, yInter, midUnc)
        if pointsInside < targetPoints:
            minUnc = midUnc
        else:
            maxUnc = midUnc
    minx = min(xvals)
    maxx = max(xvals)
    dx = maxx - minx
    # Let unc be the uncertainty for the y intercept calculated above (midUnc).
    # The ends of the regression line have y coordinates of
    # ymin = minx * slope + yInter, ymax = maxx * slope + yInter.
    # Let dy = ymax - ymin be the difference between these two.
    # The y difference between the ends of the steeper diagonal can be calculated as
    # dymax = ymax + unc - (ymin - unc) = ymax - ymin + 2 * unc = dy + 2 * unc
    # Similarly, the y difference between the shallower diagonal can be calculated as
    # dymin = ymin - unc - (ymax + unc) = ymax - ymin - 2 = unc = dy - 2 * unc
    # dymax - dymin = dy + 2 * unc - (dy - 2 * unc) = 4 * unc
    # Thus, the slope of the steeper diagonal is slopemax = dymax / dx and the slope
    # of the shallower diagonal is dymin / dx. The overall uncertainty of the slope is
    # slopeUnc = abs(slopmax - slopemin) / 2 = abs((dymax - dymin) / dx) / 2 =
    # =  abs(4 * unc / dx) / 2 = 2 * unc / dx (since both unc and dx are always positive)
    slopeUnc = 2 * midUnc / dx
    return (slopeUnc, midUnc)
