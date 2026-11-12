import os, subprocess, re, sys
from pathlib import Path

PROJECT_DIR=Path(__file__).resolve().parents[1]
MARKER=PROJECT_DIR/'tools'/'java17_path.txt'

def version(java):
    try:
        r=subprocess.run([str(java),'-version'],capture_output=True,text=True,timeout=10)
        text=(r.stderr or '')+'\n'+(r.stdout or '')
        m=re.search(r'version "([^"]+)',text)
        if not m: return None
        v=m.group(1)
        return int(v.split('.')[1]) if v.startswith('1.') else int(v.split('.')[0])
    except Exception:
        return None

def candidates():
    out=[]
    if MARKER.is_file():
        try: out.append(MARKER.read_text(encoding='utf-8').strip())
        except Exception: pass
    jh=os.environ.get('JAVA_HOME','').strip()
    if jh: out.append(os.path.join(jh,'bin','java.exe'))
    for name in ('java.exe','java'):
        try: out += subprocess.check_output(['where.exe',name],text=True,stderr=subprocess.DEVNULL).splitlines()
        except Exception: pass
    for base in [r'C:\Program Files\Eclipse Adoptium',r'C:\Program Files\Java',r'C:\Program Files\Microsoft\jdk']:
        p=Path(base)
        if p.is_dir(): out += [str(x) for x in p.rglob('java.exe')]
    seen=set()
    for c in out:
        c=os.path.abspath(os.path.expandvars(str(c)).strip('"'))
        if c in seen or not os.path.isfile(c): continue
        seen.add(c)
        v=version(c)
        if v is not None and v>=17:
            yield c,v

for java,v in candidates():
    MARKER.parent.mkdir(parents=True,exist_ok=True)
    MARKER.write_text(java,encoding='utf-8')
    print(f'Java {v}: {java}')
    sys.exit(0)
print('No Java 17+ runtime found.',file=sys.stderr)
sys.exit(1)
