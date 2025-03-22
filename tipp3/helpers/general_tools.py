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
        # implemented by Panzi:
        # https://gist.github.com/panzi/b4a51b3968f67b9ff4c99459fb9c5b3d
        lines = []
        for line_str in text.split('\n'):
            line = []
            line_len = 0
            # split line to words
            for word in line_str.split():
                word_len = len(word)
                next_len = line_len + word_len
                if line:
                    # add white space to separate with previous item
                    next_len += 1
                # wrapping text
                if next_len > width:
                    lines.append(' '.join(line))
                    line.clear()
                    line_len = 0
                elif line:
                    # add whitespace
                    line_len += 1

                line.append(word)
                line_len += word_len
            # last item
            lines.append(' '.join(line))
        return lines

        #if '\n' in text:
        #    temp = text.split('\n')
        #    ret = []
        #    for _splice in [argparse.RawDescriptionHelpFormatter._split_lines(self, x, width)
        #            for x in temp]:
        #        ret.extend(_splice)
        #    return ret
        #return argparse.RawDescriptionHelpFormatter._split_lines(self, text, width)
