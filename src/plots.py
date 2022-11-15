from matplotlib import pyplot as plt


def add_panel_text(ax, text, xplace=-0.15, fsz=17):
    ax.text(xplace, 1.1, text, transform=ax.transAxes, fontsize=fsz,
            fontweight='bold', va='top', ha='right')
