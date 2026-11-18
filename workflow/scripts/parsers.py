"""Pure parsing helpers shared by the workflow scripts; no I/O, so they are unit-tested in .tests/unit."""
import html
import os
from html.parser import HTMLParser


def db_hits(function):
    """BAGEL4 GeneTable 'function' cells glue 'score;Name' pairs together:
    '74.1;Mutacin_II(MutacinH-29B)166.2;Plantaricin_423' -> [(74.1, 'Mutacin_II(MutacinH-29B)'), (166.2, 'Plantaricin_423')].
    A numeric token scores the next token; a float glued to a name's tail belongs to the following pair."""
    def tail_float(token):
        for i in range(1, len(token)):
            rest = token[i:]
            if "." not in rest:
                continue
            try:
                return token[:i], float(rest)
            except ValueError:
                continue
        return None

    out, score = [], None
    for token in (t.strip() for t in function.split(";")):
        if score is not None:
            split = tail_float(token)
            name, nxt = split if split else (token, None)
            out.append((score, name))
            score = nxt
            continue
        try:
            score = float(token) if "." in token else None
        except ValueError:
            score = None
    return out


def acc_from_upload(name):
    """BAGEL4 drops the underscore from uploaded file names: 'GCA061326205.1.fna' -> 'GCA_061326205.1'."""
    stem = name.removesuffix(".fna")
    return "GCA_" + stem[3:] if stem.startswith("GCA") and not stem.startswith("GCA_") else stem


def accession_from_busco_json(path):
    """results/qc/busco/<accession>.json -> <accession> (the parent folder is 'busco', not the accession)."""
    return os.path.basename(path).removesuffix(".json")


def cas_type(products):
    """Cas system type from Cas protein names, by signature genes (Makarova et al. 2020)."""
    s = " ".join(products).lower()
    types = []
    if "cas3" in s or "cst" in s or "cas8b" in s or "i-b" in s:
        types.append("I-B")
    if "cas9" in s or "csn2" in s:
        types.append("II-A")
    if "csm" in s or "cas10" in s or "iii-a" in s:
        types.append("III-A")
    return " + ".join(types) or "fragment"


class _Text(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []

    def handle_data(self, data):
        self.parts.append(data)


def html_text(cell):
    """TYGS table cells carry HTML (<b>, <I>, <a>): "<b>'GCA_1.1.fna'</b> " -> 'GCA_1.1'."""
    p = _Text()
    p.feed(str(cell))
    return html.unescape("".join(p.parts)).strip().strip("'").removesuffix(".fna").strip()
