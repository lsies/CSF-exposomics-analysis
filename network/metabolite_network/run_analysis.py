from pathlib import Path
import argparse
import subprocess
import sys

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--figures-only',action='store_true')
    parser.add_argument('--online',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    scripts=['09_plot_network.py'] if args.figures_only else ['01_prepare_metabolites.py','02_download_reaction_sources.py','03_build_metabolite_network.py','04_fetch_string.py','05_retrieve_metascape.py','06_integrate_expression.py','08_summarize_network.py','09_plot_network.py','10_validate_package.py']
    for name in scripts:
        command=[sys.executable,str(root/'code'/name)]
        if args.online and name[:2] in ['02','04','05']:
            command.append('--online')
        print(name,flush=True)
        subprocess.run(command,cwd=root,check=True)

if __name__=='__main__':
    main()
