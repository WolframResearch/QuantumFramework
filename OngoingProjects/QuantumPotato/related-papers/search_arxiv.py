import urllib.request, urllib.parse, time, re, xml.etree.ElementTree as ET

QUERIES = {
 "T1 native: SIC factorization":            'all:"informationally complete" AND all:"factoriz"',
 "T2 native: SIC product distribution":     'all:"SIC-POVM" AND all:"product distribution"',
 "T3 wigner: product marginals":            'all:"discrete Wigner" AND all:marginals AND all:product',
 "T4 wigner: factorizable":                 'all:"factorizable" AND all:"Wigner function"',
 "T5 KD: uncorrelated":                     'all:"Kirkwood-Dirac" AND all:uncorrelated',
 "T6 KD: product / independence":           'all:"Kirkwood-Dirac" AND all:independence',
 "T7 stats: Matthews correlation quantum":  'cat:quant-ph AND all:"Matthews correlation"',
 "T8 abstract: independent binary marginals qubit": 'cat:quant-ph AND all:"independent" AND all:"marginals" AND all:qubit',
 "T9 their paper + follow-ups":             'all:"quantum potato"',
 "T10 probability repr: qubit coins":       'all:"probability representation" AND all:qubit AND all:coin',
}

ns = {'a': 'http://www.w3.org/2005/Atom'}
for name, q in QUERIES.items():
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode(
        {"search_query": q, "max_results": "12", "sortBy": "relevance"})
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            tree = ET.fromstring(r.read())
        print(f"\n### {name}  [{q}]")
        entries = tree.findall('a:entry', ns)
        if not entries:
            print("   (no hits)")
        for e in entries:
            eid = e.find('a:id', ns).text.split('/abs/')[-1]
            title = re.sub(r'\s+', ' ', e.find('a:title', ns).text.strip())
            year = e.find('a:published', ns).text[:4]
            print(f"   {eid}  ({year})  {title}")
    except Exception as ex:
        print(f"\n### {name}: ERROR {ex}")
    time.sleep(3)
