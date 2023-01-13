#!/usr/bin/env python
class square_subplots():
    def __init__(self, fig):
        self.fig = fig
        self.ax = self.fig.axes[0]
        self.figw,self.figh = 0,0
        self.params = [self.fig.subplotpars.left,
                       self.fig.subplotpars.right,
                       self.fig.subplotpars.top,
                       self.fig.subplotpars.bottom,
                       self.fig.subplotpars.wspace,
                       self.fig.subplotpars.hspace]
        self.rows, self.cols = self.ax.get_subplotspec().get_gridspec().get_geometry()
        self.update(None)
        self.cid = self.fig.canvas.mpl_connect('resize_event', self.update)


    def update(self, evt):
        figw,figh = self.fig.get_size_inches()
        if self.figw != figw or self.figh != figh:
            self.figw = figw; self.figh = figh
            l,r,t,b,wspace,hspace = self.params
            axw = figw*(r-l)/(self.cols+(self.cols-1)*wspace)
            axh = figh*(t-b)/(self.rows+(self.rows-1)*hspace)
            axs = min(axw,axh)
            w = (1-axs/figw*(self.cols+(self.cols-1)*wspace))/2.
            h = (1-axs/figh*(self.rows+(self.rows-1)*hspace))/2.
            self.fig.subplots_adjust(bottom=h, top=1-h, left=w, right=1-w)
            self.fig.canvas.draw_idle()
