from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text

COLORS={'teal':'#159C99','pink':'#EC629B','purple':'#A877D4','blue':'#6EA9DA','orange':'#EEA343','green':'#B2CE54','grey':'#D7DEE3','ink':'#26333D'}
FONT_SCALE=1.8
CANVAS_SCALE=1.5
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})

def save(fig, root, name):
    folder=root/'figures'
    folder.mkdir(exist_ok=True)
    (folder/name).parent.mkdir(parents=True,exist_ok=True)
    fig.set_size_inches(fig.get_size_inches()*CANVAS_SCALE,forward=True)
    for text in fig.findobj(match=Text):
        text.set_fontsize(text.get_fontsize()*FONT_SCALE)
    fig.canvas.draw()
    sizes=[text.get_fontsize() for text in fig.findobj(match=Text) if text.get_visible() and text.get_text().strip()]
    metadata_path = root/'data/figure_source/figure_typography.csv'
    records = pd.read_csv(metadata_path).to_dict('records') if metadata_path.exists() else []
    records = [r for r in records if r['figure'] != name]
    records.append({'figure':name,'font_scale':FONT_SCALE,'canvas_scale':CANVAS_SCALE,'width_inches':fig.get_figwidth(),'height_inches':fig.get_figheight(),'minimum_font_points':min(sizes),'maximum_font_points':max(sizes)})
    pd.DataFrame(records).to_csv(root/'data/figure_source/figure_typography.csv',index=False)
    for extension in ['pdf','svg','png','tiff']:
        options={'dpi':600 if extension=='tiff' else 180}
        if extension=='tiff':
            options['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(folder/f'{name}.{extension}',**options)
    plt.close(fig)
