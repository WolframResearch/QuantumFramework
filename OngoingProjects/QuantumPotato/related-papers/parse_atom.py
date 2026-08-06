import sys, re, xml.etree.ElementTree as ET
ns = {'a': 'http://www.w3.org/2005/Atom'}
tree = ET.parse(sys.argv[1])
entries = tree.getroot().findall('a:entry', ns)
if not entries:
    print("   (no hits)")
for e in entries:
    eid = e.find('a:id', ns).text.split('/abs/')[-1]
    title = re.sub(r'\s+', ' ', e.find('a:title', ns).text.strip())
    year = e.find('a:published', ns).text[:4]
    print(f"   {eid}  ({year})  {title}")
