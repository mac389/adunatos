import mygene, itertools

import numpy as np 

from awesome_print import ap 
from collections import Counter, OrderedDict
mg = mygene.MyGeneInfo()


def get_proper_gene_function_key(gene,keys,fill_value='Unknown'):
    entrezID = mg.query(gene,scopes='symbol',species='human',fields='entrezgene')['hits'][0]['entrezgene']
    try:
        ontology = mg.getgene(entrezID,fields=['go'])['go']
    except:
        return 'Unknown' #Error means protein lacks Entrez annotation
    terms = [] 
    for item in ontology.values():
        if type(item)==list:
            terms.extend([subitem['term'] for subitem in item])
        elif type(item) == dict:
            terms.append(item['term'])
    intersection = list(set(terms) & set(keys))
    if len(intersection) == 0:
        return 'Unknown'
    else: 
        return intersection.pop()
def mk_groups(data):
    try:
        newdata = data.items()
    except:
        return

    thisgroup = []
    groups = []
    for key, value in newdata:
        newgroups = mk_groups(value)
        if newgroups is None:
            thisgroup.append((key, value))
        else:
            thisgroup.append((key, len(newgroups[-1])))
            if groups:
                groups = [g + n for n, g in zip(newgroups, groups)]
            else:
                groups = newgroups
    return [thisgroup] + groups

def add_line(ax, xpos, ypos):
    line = plt.Line2D([xpos, xpos], [ypos + .1, ypos],
                      transform=ax.transAxes, color='black')
    line.set_clip_on(False)
    ax.add_line(line)

def label_group_bar(ax, data):
    groups = mk_groups(data)
    xy = groups.pop()
    x, y = zip(*xy)
    ly = len(y)
    xticks = range(1, ly + 1)

    ax.bar(xticks, y, align='center')
    ax.set_xticks(xticks)
    ax.set_xticklabels(x)
    ax.set_xlim(.5, ly + .5)
    ax.yaxis.grid(True)

    scale = 1. / ly
    for pos in xrange(ly + 1):
        add_line(ax, pos * scale, -.1)
    ypos = -.2
    while groups:
        group = groups.pop()
        pos = 0
        for label, rpos in group:
            lxpos = (pos + .5 * rpos) * scale
            ax.text(lxpos, ypos, label, ha='center', transform=ax.transAxes)
            add_line(ax, pos * scale, ypos)
            pos += rpos
        add_line(ax, pos * scale, ypos)
        ypos -= .1

def retrieve_annotation(id_list):
    request = Entrez.epost("gene",id=",".join(str(id_list)))
    try:
        result = Entrez.read(request)
        ap(result)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print "An error occurred while retrieving the annotations."
        print "The error returned was %s" % e
 
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)
 
    print "Retrieved %d annotations for %d genes" % (len(annotations),
            len(id_list))
 
    return annotations

def format(label):
    words = label.split()
    if len(words)>6 or len(label)>50:
        return ' '.join([' '.join(words[:3]),'\n',' '.join(words[3:])])
    else:
        return label

def get_gene_function(gene,levels=None):
    try:
        entrezID = mg.query(gene,scopes='symbol',species='human',fields='entrezgene')['hits'][0]['entrezgene']
        ontology = mg.getgene(entrezID,fields=['go'])['go']
        terms = [value[0]['term'] for value in ontology.values()]
        if not levels:
            return terms[-1]
        elif levels == 'all':
            return terms
        elif type(levels) == int:
            return terms
        else:
            return terms[:-levels] if levels >= len(terms) else terms
    except:
        return 'Unknown'

def get_occurences(word,seq):
  return [i for i, x in enumerate(seq) if x == word]

def word_change_boundaries(list_of_words,cutoff=8): #Still buggy
  return {word:int(np.median(get_occurences(word,list_of_words)))
          for word in set(list_of_words) if len(get_occurences(word,list_of_words))>cutoff}
  
def get_concept_boundaries(list_of_words,cutoff=8):
    current_word = list_of_words[0]
    boundaries = []
    word_frequencies = Counter(list_of_words)
    for i in xrange(1,len(list_of_words)):
        if current_word != list_of_words[i]: 
            if  word_frequencies[current_word]>cutoff:
                boundaries+=[i]
            current_word = list_of_words[i]
    return boundaries

def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key = seq.__getitem__)

def flatten(container):
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
        else:
            yield i

def array_from_lists(lst):
  nlists = len(lst)
  nwords = max([len(elem[0]) for elem in lst])
  
  #Strings have length 1
  ans = np.full((nlists,nwords),'None',dtype=object)
  for i,sublist in enumerate(lst):
    if type(sublist) == str or type(sublist) == unicode:
      ans[i,0] = sublist
    else:
      for j,word in enumerate(list(flatten(sublist))):
        ans[i,j] = word
  return ans.astype(str)

def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def uniqfy(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def pairwise(iterable): #From itertools recipies for Python 2.7
  "s -> (s0,s1), (s1,s2), (s2, s3), ..."
  a, b = itertools.tee(iterable)
  next(b, None)
  return itertools.izip(a, b)

def nearest_not(arr,idx,nullval):
  return [row[idx] if row[idx] !=nullval else row[row!=nullval][np.argmin(np.where(row!=nullval)[0]-5)]
                                              for row in arr]
  
def after(path,parent):
    return [structure for structure in path if path.index(structure) > path.index(parent)]

def before(path,parent):
    return [structure for structure in path if path.index(structure) < path.index(parent)]

def allChildrenOfParent(parent,ontology):
    paths = list(set(itertools.chain.from_iterable([after(path.split('_'),parent) 
                for path in ontology if "%s_"%parent in path])))
    return paths

def delateralize(structures_ontology):
  ans = structures_ontology[0]

  ans[-1] = ans[-1].replace(", left","").replace(", right","")
  return ans 

def get_proper_key(item,dictionary,ontology,fill_value='Unknown'):
  possible_clusters = set(dictionary.keys())
  if 'choroid plexus' in item:
    ans = fill_value #This is an ugly hack
  else:
    k = uniqfy(flatten(list(k for k,_ in itertools.groupby([path.split('_') 
      for path in ontology if item in path]))))

    if len(list(set(k) & possible_clusters)) > 0:
      return list(set(k) & possible_clusters)[0]
    else:
      return fill_value #Ugly hack, return value should reflect level of ontology in dictionary of structures
    ans = ans[0] if len(ans)>0 else item
  return ans

def ancestors(currentNode,ontology):
  ans =  list(k for k,_ in itertools.groupby([path.split('_') 
      for path in ontology if currentNode in path]))

  if len(ans) == 0:
    return currentNode
  elif len(ans) == 1:
    return list(flatten(ans))
  elif len(ans) == 2:
    return list(flatten(delateralize(ans))) #Ignoring laterality. Unify lists that differ only in laterality. 
  elif len(ans) >= 3:
    return list(flatten(ans[0])) #This approxiation is accurate to within the region

def ancestor(currentNode,ontology,levels=1): #Count of levels of ancestors to return
  paths = [path for path in ontology if "%s_"%currentNode in path] 
  return list(set(itertools.chain.from_iterable([[structure for structure in path.split('_')  
            if path.split('_').index(currentNode) > path.split('_').index(structure)] for path in paths])))[-levels:]
  #Flattening won't work structure paths of very different lengths, I don't think . 

def rcm(a,b):
  #Return lowest ontology in both
  common_ontologies = [one for one,two in itertools.product(a[::-1],b[::-1]) if one == two]
  return common_ontologies[0] if len(common_ontologies)>0 else []

def intersection(tup):
  a,b = tup
  a = uniqfy(a) if type(a) == list else uniqfy([a])
  b = uniqfy(b) if type(b) == list else uniqfy([b])

  return rcm(a,b) 