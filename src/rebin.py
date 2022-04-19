import numpy as np
import numba as nb

@nb.njit
def rebin(old_bins, old_vals, new_bins, err_check = True):
    """Rebins data defined on bins onto a new set of target bins.

    Parameters
    ----------
    old_bins : numpy.ndarray
        Edges of bins for which old_vals are defined. Length should be
        len(old_bins) - 1 = len(old_vals)
    old_vals : numpy.ndarray
        Values defined on old_bins.
    new_bins : numpy.ndarray
        Edges of target bin that you want to rebin to.
    err_check : bool
        True if you want to check inputs are valid.

    Returns
    -------
    numpy.ndarray
        Values defined on new_bins, which are the result of rebinning.
    """
    
    n_old = old_vals.shape[0]
    n_new = new_bins.shape[0] - 1
    
    if err_check:
        
        # check shape
        if n_old+1 != old_bins.shape[0]:
            raise ValueError("old_bins has the wrong length.")
        
        # check that new bins overlap with the old bins
        if new_bins[0] < old_bins[0] and new_bins[1] < old_bins[0]:
            raise ValueError("bins do not overlap.")

        if new_bins[n_new-1] > old_bins[n_old] and new_bins[n_new] > old_bins[n_old]:
            raise ValueError("bins do not overlap.")

        # check that bins are all increasing
        for i in range(n_old):
            if old_bins[i+1] <= old_bins[i]:
                raise ValueError("old_bins are not always increasing.")
                
        for i in range(n_new):
            if new_bins[i+1] <= new_bins[i]:
                raise ValueError("new_bins are not always increasing.")
          
    # Do binning
    new_vals = np.zeros(n_new)
    l = 0
    
    for i in range(n_new):
        b1 = new_bins[i+1] - new_bins[i]
        for j in range(l,n_old):
            # four different cases
        
            # ______       (old bin)
            #     ________ (new bin)
            if old_bins[j] <= new_bins[i] and \
               old_bins[j+1] > new_bins[i] and old_bins[j+1] < new_bins[i+1]:
          
                b2 = old_bins[j+1] - new_bins[i]
                new_vals[i] = new_vals[i] + (b2/b1)*old_vals[j]
          
            #    ____    (old bin)
            #  ________  (new bin)
            elif old_bins[j] >= new_bins[i] and old_bins[j+1] <= new_bins[i+1]:
          
                b2 = old_bins[j+1] - old_bins[j]
                new_vals[i] = new_vals[i] + (b2/b1)*old_vals[j]
          
            #        ____    (old bin)
            #  ________      (new bin)
            elif old_bins[j] >= new_bins[i] and \
                 old_bins[j] < new_bins[i+1] and old_bins[j+1] > new_bins[i+1]:
        
                b2 = new_bins[i+1] - old_bins[j]
                new_vals[i] = new_vals[i] + (b2/b1)*old_vals[j]
        
            # ________    (old bin)
            #   ____      (new bin) 
            elif old_bins[j] < new_bins[i] and old_bins[j+1] > new_bins[i+1]:
          
                new_vals[i] = new_vals[i] + old_vals[j]
          
            #         _____   (old bin)
            #   ____          (new bin) 
            elif old_bins[j] > new_bins[i+1]:
                # time to move onto the next new_bin
                l = j - 1
                break
            else:
                pass
                # nothing
    
    return new_vals