import re
from enum import Enum

from app.tools.logger import clinical_trials_logger

logger = clinical_trials_logger()


class BulletPoint(Enum):
    ORDERED = "."
    STAR = "*"
    ORDERED_BRACKETS = ")"
    NONE = ""


class Line:
    ordered_bracket_pattern = r"^\[(\d+)\]"
    star_pattern = r"^\*"
    ordered_pattern = r"^\d+\."

    def __init__(self, line: str):
        logger.debug(f"Initialising Line class for {line}")
        self.leading_spaces: int = Line.count_leading_spaces_and_tabs(line)
        self.line: str = line.strip()
        self.type: BulletPoint = Line.type(self.line)
        logger.debug(f"leading_space = {self.leading_spaces}, line = {self.line}, type = {self.type}")

    def __str__(self):
        return (f"Line: leading_spaces={self.leading_spaces}, "
                f"line={self.line}, type={self.type}")

    def __repr__(self):
        return f"Line({self.leading_spaces!r}, {self.line!r}, {self.type!r})"

    @staticmethod
    def count_leading_spaces_and_tabs(line: str):
        leading_whitespace = len(line) - len(line.lstrip(' \t'))
        return leading_whitespace

    @staticmethod
    def type(new_line: str):
        if re.match(Line.star_pattern, new_line) is not None:
            return BulletPoint.STAR
        if re.match(Line.ordered_bracket_pattern, new_line) is not None:
            return BulletPoint.ORDERED_BRACKETS
        if re.match(Line.ordered_pattern, new_line) is not None:
            return BulletPoint.ORDERED
        return BulletPoint.NONE

    def number(self):
        if self.type == BulletPoint.NONE:
            return -1
        if self.type == BulletPoint.STAR:
            return 0
        if self.type == BulletPoint.ORDERED:
            return int(self.line.split(".")[0])
        if self.type == BulletPoint.ORDERED_BRACKETS:
            return int(self.line.split(")")[0][1:])

    def text(self):
        if self.type == BulletPoint.NONE:
            return self.line
        start_index = self.line.find(self.type.value)
        if start_index == -1:
            raise ValueError("stuff")
        if start_index == len(self.line) - 1:
            raise ValueError("different stuff")
        return self.line[start_index + 1:].strip()

    def html(self):
        if self.type == BulletPoint.NONE:
            return "<p>{text}</p>".format(text=self.text())
        return "<li>{text}</li>".format(text=self.text())


def list_start_finish(line_type: BulletPoint, start=True):
    logger.debug(f"eligibility.list_starter_finisher({line_type}, {start})")
    if line_type == BulletPoint.NONE:
        raise ValueError("If line type is NONE, there is no start or finish because it is not a list")
    if start:
        if line_type == BulletPoint.STAR:
            logger.debug("html needed is <ul>")
            return "<ul>"
        else:
            logger.debug("html needed is <ol>")
            return "<ol>"
    if line_type == BulletPoint.STAR:
        logger.debug("html needed is </ul>")
        return "</ul>"
    else:
        logger.debug("html needed is </ol>")
        return "</ol>"


def is_same(old_line: Line, new_line: Line):
    logger.debug(f"elibility.is_same({old_line}, {new_line})")
    if old_line.type != new_line.type:
        return False
    if old_line.leading_spaces != new_line.leading_spaces:
        return False
    if new_line.type == BulletPoint.ORDERED or new_line.type == BulletPoint.ORDERED_BRACKETS:
        if new_line.number == 1:
            return False
    logger.debug("The lines are at the same bullet point level")
    return True


def format_eligibility(eligibility_string: str):
    logger.debug(f"eligibility.format_eligibility({eligibility_string})")
    lines = eligibility_string.splitlines()
    html_list = []
    type_nest = []  # list of all the types of bullet points currently nested
    number_nest = []  # list of all the spaces of bullet point sections currently nested
    old_line = None
    for i in range(0, len(lines)):
        new = Line(lines[i])
        if new.line == "":
            continue
        if old_line:  # i.e. isn't the start of a new top level bullet point or a random sentence
            if not is_same(old_line, new) or new.number() == 1 or old_line.leading_spaces != new.leading_spaces:
                if new.type == BulletPoint.NONE:
                    logger.debug("All bullet point blocks have ended")
                    html_list += [list_start_finish(line_type, start=False) for line_type in list(reversed(type_nest))]
                    type_nest, number_nest = [], []
                elif (new.leading_spaces > old_line.leading_spaces or new.number() == 1 or
                      (new.leading_spaces == old_line.leading_spaces and new.number() == 0)):
                    logger.debug("Starting a new bullet point block that's further indented from the current")
                    html_list.append(list_start_finish(new.type, start=True))
                    number_nest.append(new.leading_spaces)
                    type_nest.append(new.type)
                else:  # going up
                    for ind in range(len(type_nest) - 1, -1, -1):
                        logger.debug("Calculating which level of indent we're retreating to")
                        line_type = type_nest[ind]
                        space = number_nest[ind]
                        if line_type != new.type or space > new.leading_spaces:
                            logger.debug("Retreating a level")
                            html_list.append(list_start_finish(line_type, start=False))
                            type_nest.pop()
                            number_nest.pop()
                        if space < new.leading_spaces:  # seems unlikely
                            # FIXME what is the correct response if we fail to parse an eligibility string?
                            raise ValueError("Bullet Points are not parseable")
                        else:
                            # space == new.leading_spaces, type == new.type
                            logger.debug("Reached correct level")
                            break
        else:
            if new.type != BulletPoint.NONE:
                logger.debug("Starting a new bullet point block from scratch")
                html_list.append(list_start_finish(new.type, start=True))
        logger.debug("Finished calculating bullet point positioning, adding line to html list")
        html_list.append(new.html())
        old_line = new
    logger.debug("Retreating out of all levels of indent")
    html_list += [list_start_finish(i, False) for i in list(reversed(type_nest))]
    return html_list
