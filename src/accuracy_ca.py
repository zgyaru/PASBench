import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn import datasets
from sklearn.model_selection import KFold,StratifiedKFold
from collections import  Counter
from sklearn.metrics import accuracy_score
import pandas as pd
import seaborn as sns
from itertools import product
from sklearn.metrics import confusion_matrix
from sklearn.utils import check_matplotlib_support
#from sklearn.utils.validation import _deprecate_positional_args
from sklearn.base import is_classifier
from sklearn import preprocessing
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
import warnings
import os,glob
from  scipy.stats import ttest_rel
from scipy.stats.mstats import gmean
warnings.filterwarnings('ignore')
sns.set_context("paper",font_scale=1.8)


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    #cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    #cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="left",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=1, markeredgecolor ='b',fillstyle='full')
    ax.tick_params(which="minor", bottom=False, left=False, right=False)

    return im


def annotate_heatmap(im, data=None, dataType = 'time',
                     valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, NA_text = 'time', **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = np.nanmin(data)+(np.nanmax(data)-np.nanmin(data))/2

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    if dataType == 'time':
        texts = []
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                tt = 0
                if np.isnan(data[i,j]):
                    text = im.axes.text(j, i, " - ", **kw)
                    texts.append(text)
                else:
                    if data[i, j]>60:
                        if data[i, j]>3600:
                            tt = data[i, j]/60/60
                            valfmt="{x:.1f} h"
                        else:
                            tt = data[i, j]/60
                            valfmt="{x:.1f} m"      
                    else:
                        tt = data[i, j]
                        valfmt="{x:.1f} s"
                        
                    if isinstance(valfmt, str):
                        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

                    kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                    text = im.axes.text(j, i, valfmt(tt), **kw)
                    texts.append(text)
    elif dataType == 'acc':
        texts = []
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                tt = 0
                if np.isnan(data[i,j]):
                    kw.update(color='r')
                    text = im.axes.text(j, i, NA_text, **kw)
                    texts.append(text)
                else:
                    tt = data[i, j]
                    valfmt="{x:.2f}"
                
                    if isinstance(valfmt, str):
                        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

                    kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                    text = im.axes.text(j, i, valfmt(tt), **kw)
                    texts.append(text)
    
    elif dataType == 'memory':
        texts = []
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                tt = 0
                if np.isnan(data[i,j]):
                    text = im.axes.text(j, i, " - ", **kw)
                    texts.append(text)
                else:
                    if data[i, j]>1024:
                        tt = data[i, j]/1024
                        valfmt="{x:.1f} G"      
                    else:
                        tt = data[i, j]
                        valfmt="{x:.1f} M"
                        
                    if isinstance(valfmt, str):
                        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

                    kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                    text = im.axes.text(j, i, valfmt(tt), **kw)
                    texts.append(text)

    return texts

def parseY(Y):
    dic = {}
    i=0
    for y in Y:
        if y in dic:
            pass
        else:
            dic[y] = i
            i += 1
    rr = []
    for y in Y:
        rr.append(dic[y])
    return(np.array(rr))

def get_kfold(X,Y,k=5,model='lr',l=1e3, scale='minMax',n_estimators=100):
    ## scale: 'none', 'minMax', 'stand'
    y_true_all = np.array([])
    y_pred_all = np.array([])
    kf = StratifiedKFold(n_splits=k,shuffle=True)
    acc = []
    for train_index, test_index in kf.split(X, Y):
        x_train = X[train_index]
        x_test = X[test_index]
        y_train = Y[train_index]
        y_test = Y[test_index]
        
        if scale == 'minMax':
            scaler = preprocessing.MinMaxScaler().fit(x_train)
            x_train = scaler.transform(x_train)
            x_test = scaler.transform(x_test)
        elif scale == 'stand':
            scaler = preprocessing.StandardScaler().fit(x_train)
            x_train = scaler.transform(x_train)
            x_test = scaler.transform(x_test)
        elif scale == 'none':
            pass
        
        if model == 'lr':
            clr = LogisticRegression(C=l)
        elif model == 'rfc':
            clr = RandomForestClassifier(n_estimators)
        elif model == 'svc':
            clr = LinearSVC(C= l)
        clr.fit(x_train, y_train)
        acc.append(clr.score(x_test,y_test))
        #y_true_all = np.concatenate([y_true_all, y_test])
        #y_pred_all = np.concatenate([y_pred_all, clr.predict(x_test)])
    return(np.mean(acc))





