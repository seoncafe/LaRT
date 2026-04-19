    def slice_plot(self, axis='z', value=None, quantity='nH',
                   ax=None, log=True, cmap='viridis',
                   show_leaf_boundaries=False,
                   boundary_color='k', boundary_lw=0.3,
                   boundary_alpha=0.7,
                   show_leaf_centers=False,
                   center_color='k', center_marker='o',
                   center_size=8, center_alpha=0.9, **kwargs):
        """
        Quick 2-D slice plot of a physical quantity through the grid.

        Parameters
        ----------
        axis : {'x', 'y', 'z'}
            Normal axis of the slice plane.
        value : float, optional
            Position of the slice along ``axis``.  Defaults to box centre.
        quantity : {'nH', 'T', 'vx', 'vy', 'vz'}
            Quantity to plot.
        ax : matplotlib.axes.Axes, optional
        log : bool
            Plot log10 of the quantity if True.
        cmap : str
        show_leaf_boundaries : bool
            If True, overlay boundaries of the leaf cells that intersect
            the slice plane.
        boundary_color : str
            Edge colour for the leaf-cell boundary overlay.
        boundary_lw : float
            Line width for the leaf-cell boundary overlay.
        boundary_alpha : float
            Transparency for the leaf-cell boundary overlay.
        show_leaf_centers : bool
            If True, overlay markers at the centres of leaf cells that
            intersect the slice plane.
        center_color : str
            Marker colour for the leaf-cell centre overlay.
        center_marker : str
            Matplotlib marker style for the leaf-cell centre overlay.
        center_size : float
            Marker size for the leaf-cell centre overlay.
        center_alpha : float
            Transparency for the leaf-cell centre overlay.
        **kwargs
            Passed to ``matplotlib.patches.Rectangle`` for the filled cells.

        Returns
        -------
        im : matplotlib.collections.PatchCollection
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.collections import PatchCollection
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        if ax is None:
            _, ax = plt.subplots()

        if value is None:
            value = self.origin[{'x': 0, 'y': 1, 'z': 2}[axis]] + self.boxlen * 0.5

        ax_map = {'x': ('y', 'z', 'cy', 'cz', 'cx'),
                  'y': ('x', 'z', 'cx', 'cz', 'cy'),
                  'z': ('x', 'y', 'cx', 'cy', 'cz')}
        ha, va, ca, da, na = ax_map[axis.lower()]

        patches, values = [], []
        boundary_patches = []
        center_x, center_y = [], []
        for lf in self.leaves():
            c_normal = getattr(lf, na)
            if abs(c_normal - value) <= lf.h:
                cx_ = getattr(lf, ca)
                cy_ = getattr(lf, da)
                rect = mpatches.Rectangle(
                    (cx_ - lf.h, cy_ - lf.h), 2*lf.h, 2*lf.h, **kwargs)
                patches.append(rect)
                if show_leaf_boundaries:
                    boundary_patches.append(
                        mpatches.Rectangle(
                            (cx_ - lf.h, cy_ - lf.h), 2*lf.h, 2*lf.h
                        )
                    )
                if show_leaf_centers:
                    center_x.append(cx_)
                    center_y.append(cy_)
                v = getattr(lf, quantity)
                values.append(np.log10(max(v, 1e-100)) if log else v)

        pc = PatchCollection(patches, cmap=cmap, linewidths=0)
        pc.set_array(np.array(values))
        ax.add_collection(pc)

        if show_leaf_boundaries and boundary_patches:
            pc_boundary = PatchCollection(
                boundary_patches,
                facecolor='none',
                edgecolor=boundary_color,
                linewidths=boundary_lw,
                alpha=boundary_alpha,
            )
            ax.add_collection(pc_boundary)

        if show_leaf_centers and center_x:
            ax.scatter(
                center_x, center_y,
                c=center_color, marker=center_marker, s=center_size,
                alpha=center_alpha
            )

        ax.set_xlim(self.origin[0], self.origin[0] + self.boxlen)
        ax.set_ylim(self.origin[1], self.origin[1] + self.boxlen)
        ax.set_aspect('equal')
        ax.set_xlabel(ha)
        ax.set_ylabel(va)
        label = f'log10({quantity})' if log else quantity

        cax = inset_axes(
            ax,
            width="4%",
            height="100%",
            loc='lower left',
            bbox_to_anchor=(1.02, 0.0, 1.0, 1.0),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
        ax.figure.colorbar(pc, cax=cax, label=label)
        return pc
