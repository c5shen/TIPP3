import psutil, os
import argparse

tqdm_styles = {
        'desc': '\tRunning...', 'ascii': False,
        'ncols': 80,
        #'disable': True,
        #'colour': 'green',
        'mininterval': 0.5
        }

# return memory usage of python process by MB
def memoryUsage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)
    return mem

# reformat argparse help text formatting
class SmartHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def add_text(self, text):
        if text is not None:
            text = text.replace("\\n", "\n").replace("\\t", "\t")
        super().add_text(text)
    def _split_lines(self, text, width):
        if '\n' in text:
            temp = text.split('\n')
            ret = []
            for _splice in [argparse.RawDescriptionHelpFormatter._split_lines(self, x, width)
                    for x in temp]:
                ret.extend(_splice)
            return ret
        return argparse.RawDescriptionHelpFormatter._split_lines(self, text, width)
