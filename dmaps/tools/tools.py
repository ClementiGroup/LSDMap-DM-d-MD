import random
import numpy as np
from scipy.sparse import spdiags

def smooth2a(arrayin, nr, nc):

    # Building matrices that will compute running sums.  The left-matrix, eL,
    # smooths along the rows.  The right-matrix, eR, smooths along the
    # columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
    # (2*Nc+1) rectangle centered on element "i".


    row = arrayin.shape[0]
    col = arrayin.shape[1]

    el = spdiags(np.ones((2*nr+1, row)),range(-nr,nr+1), row, row).todense()
    er = spdiags(np.ones((2*nc+1, col)), range(-nc,nc+1), col, col).todense()

    # Setting all "nan" elements of "arrayin" to zero so that these will not
    # affect the summation.  (If this isn't done, any sum that includes a nan
    # will also become nan.)

    a = np.isnan(arrayin)
    arrayin[a] = 0.

    # For each element, we have to count how many non-nan elements went into
    # the sums.  This is so we can divide by that number to get a mean.  We use
    # the same matrices to do this (ie, "el" and "er").

    nrmlize = el.dot((~a).dot(er))
    nrmlize[a] = None

    # Actually taking the mean.

    arrayout = el.dot(arrayin.dot(er))
    arrayout = arrayout/nrmlize

    return arrayout


def do_hist2D(x, y, nbins, nextrabins=0):

    # build a 2D histogram
    min_x = np.amin(x)
    max_x = np.amax(x)

    min_y = np.amin(y)
    max_y = np.amax(y)

    step_x = (max_x - min_x)/nbins
    step_y = (max_y - min_y)/nbins

    bins1 = min_x - (nextrabins*step_x) + np.array(range(nbins+2*nextrabins+1)) * step_x
    bins2 = min_y - (nextrabins*step_y) + np.array(range(nbins+2*nextrabins+1)) * step_y

    # rescale the last non-extras bin to include the very last point
    bins1[-1-nextrabins] += 1e-6
    bins2[-1-nextrabins] += 1e-6

    inds1 = np.digitize(x, bins1) - 1  # create array with the number of the bin every point is in
    inds2 = np.digitize(y, bins2) - 1

    hist_idxs = [[[] for idx in range(nbins+2*nextrabins)] for jdx in range(nbins+2*nextrabins)]

    # rescale back the last non-extra bins
    bins1[-1-nextrabins] -= 1e-6
    bins2[-1-nextrabins] -= 1e-6

    # create histogram
    for idx, [ind1, ind2] in enumerate(np.array([inds1, inds2]).T):
        hist_idxs[ind1][ind2].append(idx)

    return bins1, bins2, hist_idxs


def draw_points_hist2D(hist_idxs, nbins, npoints, border_frac=0.0):

    drawn_points = []

    hist = np.zeros((nbins,nbins), dtype ='int')
    for idx in range(nbins):
        for jdx in range(nbins):
            hist[idx][jdx] = len(hist_idxs[idx][jdx])

    if border_frac > 0.0:
        nborder = int(border_frac*npoints)
    
        #finding border bins
        border_bins = []
        tot_points_border = 0
        for idx in range(nbins):
            # looking from left 
            for jdx in range(nbins):
                if hist[idx][jdx] > 0:
                    which_bin = [idx, jdx]
                    if which_bin not in border_bins:
                        border_bins.append(which_bin)
                        tot_points_border += hist[which_bin[0], which_bin[1]]
                    break
            for jdx in range(nbins-1,-1,-1):
	        if hist[idx][jdx] > 0:
                    which_bin = [idx, jdx]
                    if which_bin not in border_bins:
                        border_bins.append(which_bin)
                        tot_points_border += hist[which_bin[0], which_bin[1]]
                    break
            for jdx in range(nbins):
                if hist[jdx][idx] > 0:
                    which_bin = [jdx, idx]
                    if which_bin not in border_bins:
                        border_bins.append(which_bin)
                        tot_points_border += hist[which_bin[0], which_bin[1]]
                    break
            for jdx in range(nbins-1,-1,-1):
	        if hist[jdx][idx] > 0:
                    which_bin = [jdx, idx]
                    if which_bin not in border_bins:
                        border_bins.append(which_bin)
                        tot_points_border += hist[which_bin[0], which_bin[1]]
                    break

        # draw points
        # if points on border are less than desired number draw them all
        if tot_points_border <= nborder:
            for which_bin in border_bins:
                for point in hist_idxs[which_bin[0]][which_bin[1]]:
                    drawn_points.append(point)
            npicked = tot_points_border
        else:
            npicked = 0
            while npicked < nborder:
                # draw a bin and a point in it
                which_bin = random.choice(border_bins)
                point = random.choice(hist_idxs[which_bin[0]][which_bin[1]])
    
                # if point has already been drawn take another one in the same bin but check not all the points in the bin are already in drawn_points. Next two commented lines don't work, need to fix them
                #while (point in drawn_points) and (hist_idxs[which_bin[0]][which_bin[1]] not in drawn_points):
                #    point = rd.choice(hist_idxs[which_bin[0]][which_bin[1]])
                if point not in drawn_points:
                    drawn_points.append(point)
                    npicked += 1
    else:
        npicked = 0

    # construct list with nonempty bins to save time when drawing if many bins are empty (useful?)
    nonempty_bins = []
    for idx in range(nbins):
        for jdx in range(nbins):
            if hist[idx][jdx] != 0:
                nonempty_bins.append([idx, jdx])
    # draw remaining points uniformly
    while npicked < npoints:
        which_bin = random.choice(nonempty_bins)
        point = random.choice(hist_idxs[which_bin[0]][which_bin[1]])
        if point not in drawn_points:
            drawn_points.append(point)
            npicked += 1

    return drawn_points

### does not handle empty bins at the borders
#def draw_points_hist2D_test(hist, nbins, npoints, sampling='uniform', border_frac=0.0):
#
#    drawn_points = []
#
#    if border_frac > 0.0:
#        nborder = int(border_frac*npoints)
#        borders = hist[0] + hist[-1] + [hist[idx][0] for idx in range(1,nbins-1)] + [hist[idx][-1] for idx in range(1,nbins-1)]
#        # keep only non empty bins
#        borders = [bin for bin in borders if bin]
#
#        npoints_borders = 0
#        for bin in borders:
#            npoints_borders += len(bin)
#
#        print npoints_borders, nborder
#
#        if npoints_borders <= nborder:
#            drawn_points = [point for bin in borders for point in bin]
#            npicked = npoints_borders
#        else:
#            npicked = 0
#            drawn_points = []
#            while npicked < nborder:
#                point = random.choice(random.choice(borders))
#                if point not in drawn_points:
#                    drawn_points.append(point)
#                    npicked = npicked + 1
#    else:
#        nborder = 0
#        npicked = 0
#
#    # construct list with nonempty bins to save time when drawing if many bins are empty
#    nonempty_bins = [bin for line in hist for bin in line if bin]
#
#    # draw remaining points uniformly
#    while npicked < npoints:
#        point = random.choice(random.choice(nonempty_bins))
#        if point not in drawn_points:
#            drawn_points.append(point)
#            npicked = npicked + 1
#
#    return drawn_points


#def draw_idxs(hist, nbins):
#   ind1 = random.randrange(nbins)
#   ind2 = random.randrange(nbins)
#
#   while not hist[ind1][ind2]:
#       ind1 = random.randrange(nbins)
#       ind2 = random.randrange(nbins)
#
#   return ind1, ind2
