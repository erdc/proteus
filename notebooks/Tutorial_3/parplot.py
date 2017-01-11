"""parallel plotting with IPython.parallel"""

import sys
from io import BytesIO

import matplotlib
import numpy as np
import PIL.Image

from IPython import get_ipython
from IPython.display import display, clear_output, Image


# pil2png, png2array from https://github.com/minrk/ipython_extensions
# used under BSD license

def pil2png(img):
    """convert a PIL Image to png bytes"""
    fp = BytesIO()
    img.save(fp, format='PNG')
    return fp.getvalue()


def png2pil(png):
    """Convert png bytes to a PIL Image object"""
    return PIL.Image.open(BytesIO(png))


def png2array(png):
    """Convert png bytes to a numpy uint8 RGBA array"""
    img = png2pil(png)
    arr = np.fromstring(img.tobytes(), dtype=np.uint8)
    return arr.reshape(img.size[1], img.size[0], 4)


def stack_images(pngs):
    """Stack one or more images with PIL
    
    Uses Image.alpha_composite to layer images on top of each other.
    
    Takes one arg: an iterable of png data
    
    Returns: png data of the overlayed image
    """
    
    it = iter(pngs)
    img = png2pil(next(it))
    for overlay in map(png2pil, it):
        img = PIL.Image.alpha_composite(img, overlay)
    return pil2png(img)


def _plot_f(method, *args, **kwargs):
    """function to be applied by parallel_plot"""
    if 'matplotlib.backends' not in sys.modules:
        matplotlib.use('agg')
    import matplotlib.pyplot as plt

    # make a blank figure with no margins
    fig = plt.figure(frameon=False)
    fig.subplotpars.bottom = 0
    fig.subplotpars.top = 1
    fig.subplotpars.left = 0
    fig.subplotpars.right = 1
    ax = fig.gca()
    ax.axis('off')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    # do the plot
    xlim = kwargs.pop('xlim', None)
    ylim = kwargs.pop('ylim', None)
    plot = getattr(ax, method)
    plot(*args, **kwargs)
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    
    # return the figure as transparent png
    buf = BytesIO()
    fig.canvas.print_figure(buf, transparent=True, format='png')
    plt.close(fig)
    return buf.getvalue()


def parallel_plot(view, method, *args, **kwargs):
    """Compute a matplotlib plot in parallel
    
    This calls multiple matplotlib plot commands in parallel,
    gathers the resulting PNG data, and compiles it together.
    
    This will only work properly if the plot calls leave
    the areas they are not plotting as transparent.
    
    Parameters
    ----------
    
    view: an IPython.parallel DirectView
    method: str
        The plot method name, e.g. 'plot' or 'scatter'
    {x,y}lim: (int, int), kwarg-only
        Specify the x and y limits.
        xlim, ylim typically must be specified for the split 
        This typically must be specified to ensure consistent
    imshow: bool, kwarg-only
        imshow=True uses matplotlib imshow to display the plot.
        This will display the proper axes and ticks,
        at the expense of sometimes unattractive interpolation.
    
    Any further args, kwargs are passed unmodified to the remote plot calls.
    
    The gathered result is stitched together and displayed with imshow.
    """
    
    imshow = kwargs.pop('imshow', False)
    ar = view.apply_async(_plot_f, method, *args, **kwargs)
    
    if imshow:
        import matplotlib.pyplot as plt
        fig = plt.gcf()
        ax = plt.gca()
        xlim = list(kwargs.get('xlim', []))
        ylim = list(kwargs.get('ylim', []))
        if xlim and ylim:
            extent = xlim + ylim
        else:
            extent = None
    
    pngiter = iter(ar)
    merged = next(pngiter)
    
    # show the first frame
    if imshow:
        img = png2pil(merged)
        aspect = float(img.size[1]) / img.size[0]
        result = ax.imshow(png2array(merged), extent=extent, aspect=aspect)
    else:
        display(Image(data=merged, format='png'))
    
    # overlay and draw each subsequent frame
    for frame in pngiter:
        merged = stack_images([merged, frame])
        clear_output(wait=True)
        if imshow:
            result = ax.imshow(png2array(merged), extent=extent, aspect=aspect)
            display(fig)
        else:
            display(Image(data=merged, format='png'))
    
    if imshow:
        plt.close(fig)
        return result
    else:
        return merged

