from __future__ import annotations

import re
import sys
from collections.abc import Iterator


TAB = "    "


def c_doc_to_rst(fpath_or_string: str, modulename: str) -> str:
    try:
        with open(fpath_or_string) as f:
            source = f.read()
    except FileNotFoundError:
        source = fpath_or_string

    rst_blocks = []
    indent = 0
    for block in iter_docstring(source):
        rst_lines = format_doc_block(block)
        rst_block = assemble_doc_block(rst_lines, indent).replace(
            f".. py:class::\n\n{TAB}",
            f".. py:class:: {modulename}."
        )
        indent = 1
        rst_blocks.append(rst_block.replace(r"\n", ""))
    return "\n\n\n".join(rst_blocks)


def iter_docstring(source: str) -> Iterator[str]:
    while True:
        idx_offset = source.find("PyDoc_STRVAR")
        source = source[idx_offset:]
        idx_start = source.find("(") + 1
        if idx_start == 0:
            return
        idx_end = source.find(");")

        yield source[idx_start:idx_end]
        source = source[idx_end:]


def remove_whitespace(source: str) -> str:
    cleaned = source.replace('"', "").strip()
    # remove any repeated new lines
    cleaned = re.sub("\n\n+", "\n\n", cleaned)
    cleaned = cleaned.replace(r"\n\n", "\n" + r"\n")
    cleaned = cleaned.strip("\n")
    # remove any indentation
    return "\n".join(line.lstrip(" \t") for line in cleaned.split("\n"))


def get_head_body(source: str) -> tuple[str, str]:
    head, body = source.split("--")
    return head.strip("\n"), body.strip("\n")


def format_head(head: str) -> list[str]:
    directive, *lines = head.split("\n")

    directive = directive.strip("/ ")
    options = []
    while lines[0].startswith("//"):
        line = lines.pop(0)
        options.append(line.strip("/ "))

    signature = lines[-1].strip(r"\n")
    signature = clean_signature(signature)

    rst_lines = [f"{directive} {signature}"]
    for option in options:
        rst_lines.append(f"{TAB}{option}")
    return rst_lines


def clean_signature(signature: str) -> str:
    name, args = signature.split("(")
    args, ret = args.split(")")
    args = [arg.strip() for arg in args.split(",")]
    arg_str = ", ".join(args)
    return f"{name}({arg_str}){ret}"


def format_doc_block(block: str) -> str:
    cleaned = remove_whitespace(block)
    try:
        head, body = get_head_body(cleaned)
        rst_lines = format_head(head)
    except ValueError:
        rst_lines = [line.strip("/ ") for line in cleaned.split("\n") if "// " in line]
        rst_lines = [rst_lines[0], *(f"{TAB}{line}" for line in rst_lines[1:]), ""]
        body = [line for line in cleaned.split("\n") if "//" not in line]
        body = "\n".join(body[1:])
    rst_lines.extend(f"{TAB}{line}" for line in body.split("\n"))
    return rst_lines


def assemble_doc_block(rst_lines: list[str], indent: int):
    tab = indent * TAB
    return tab + f"\n{tab}".join(rst_lines)
