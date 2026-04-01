import argparse


class SmartHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom help formatter that handles embedded newlines and tabs."""

    def add_text(self, text):
        if text is not None:
            text = text.replace("\\n", "\n").replace("\\t", "\t")
        super().add_text(text)

    def _split_lines(self, text, width):
        lines = []
        for line_str in text.split('\n'):
            line = []
            line_len = 0
            for word in line_str.split():
                word_len = len(word)
                next_len = line_len + word_len
                if line:
                    next_len += 1
                if next_len > width:
                    lines.append(' '.join(line))
                    line.clear()
                    line_len = 0
                elif line:
                    line_len += 1
                line.append(word)
                line_len += word_len
            lines.append(' '.join(line))
        return lines
