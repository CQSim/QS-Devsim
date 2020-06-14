# Copyright 2016 Devsim LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from devsim import *
import sys
import math


def print_header(fh):
    '''
        Write header for backround mesh view
    '''
    fh.write('View "background mesh" {\n')

def print_footer(fh):
    '''
        Write footer for backround mesh view
    '''
    fh.write('};\n')

def get_edge_index(device, region):
    '''
    maps element edges to regular edges
    '''
    # now iterate over the edges of the element
    element_model(device=device, region=region,
                  name="eindex", equation="edge_index")
    eindex = get_element_model_values(
        device=device, region=region, name='eindex')
    eindex = [int(x) for x in eindex]
    return eindex

def get_node_index(device, region):
    '''
    maps head and tail nodes of from their edge index
    '''
    # identify all edges that need to be bisected
    # ultimately translated to an element
    edge_from_node_model(node_model="node_index", device=device, region=region)
    nindex = list(
      zip(
        [int(x) for x in get_edge_model_values(device=device, region=region, name="node_index@n0")],
        [int(x) for x in get_edge_model_values(device=device, region=region, name="node_index@n1")],
      )
    )
    return nindex

def calculate_clengths(device, region, model_values):
    '''
    calculate the characteristic lengths for each edge by bisecting the edge length
    '''
    clengths = list(get_edge_model_values(device=device, region=region, name="EdgeLength"))
    bisection_count = 0
    for i, v in enumerate(model_values):
        if v != 0:
            clengths[i] *= 0.5
            bisection_count += 1
    print("Edge Bisections: %d" % bisection_count)
    return clengths


def get_output_elements3(nindex, eindex, clengths, number_nodes, mincl, maxcl):
    '''
    gets the node indexes and the characterisic lengths for each element
    nindex : from get_node_index
    eindex : from get_edge_index
    clengths : from calculate_clengths
    number_nodes : number of nodes
    mincl : minimum characteristic length
    maxcl : maximum characteristic length
    '''
    # set upper limit to maxcl
    node_map = [maxcl] * number_nodes
    # get node indexes for each edge
    for i, n in enumerate(nindex):
        # clip minimum value to mincl
        v = max(clengths[i], mincl)
        for ni in n:
            node_map[ni] = min(node_map[ni], v) 

    #break into a per element basis
    outputelements = []
    for i in range(0, len(eindex), 3):
        ndict = {}
        # mapping of element edge into an edge index
        for j in eindex[i:i+3]:
            # mapping of edge index into a node index
            for k in nindex[j]:
                if k not in ndict:
                    ndict[k] = node_map[k]
        outputelements.append(tuple(ndict.items()))
    return outputelements

def print_elements(fh, device, region, elements):
    '''
    print background mesh triangles
    '''
    x = get_node_model_values(device=device, region=region, name="x")
    y = get_node_model_values(device=device, region=region, name="y")

    for e in elements:
      coords = []
      values = []
      for n, v in e:
        coords.extend((x[n], y[n], 0.0))
        values.append(v)
      coordstring = ", ".join([format(x, "1.15g") for x in coords])
      valuestring = ", ".join([format(x, "1.15g") for x in values])
      fh.write("ST(%s) {%s};\n" % (coordstring, valuestring))

def refine_common(fh, device, region, model_values, mincl, maxcl):
    '''
    prints out the refined elements
    model_values : non-zero for edges to be bisected
    mincl : minimum characteristic length
    maxcl : maximum characteristic length
    '''
    clengths = calculate_clengths(device=device, region=region, model_values=model_values)
  
    eindex = get_edge_index(device, region)
    nindex = get_node_index(device, region)
    number_nodes = len(get_node_model_values(device=device, region=region, name="node_index"))


    outputelements = get_output_elements3(nindex=nindex, eindex=eindex, clengths=clengths, number_nodes=number_nodes, mincl=mincl, maxcl=maxcl)
    print_elements(fh=fh, device=device, region=region, elements=outputelements)

def get_oxide_model_values(device, region):
    '''
    returns a model for non-refinement
    mincl : minimum characteristic length
    maxcl : maximum characteristic length
    '''
    test_model = [0.0] * len(get_edge_model_values(device=device, region=region, name="EdgeLength"))
    return test_model

def refine_oxide_region(fh, device, region, mincl, maxcl):
    '''
    refinement for oxide regions
    mincl : minimum characteristic length
    maxcl : maximum characteristic length
    '''
    # apply no refinement
    test_model = get_oxide_model_values(device=device, region=region)
    refine_common(fh=fh, device=device, region=region, model_values=test_model, mincl=mincl, maxcl=maxcl)

def get_silicon_model_values(device, region):
    '''
    returns a model for refinement of silicon regions
    '''
    # edge to node mapping (node0, node1)
    node_index = get_node_index(device=device, region=region)

    potential = get_node_model_values(device=device, region=region, name="Potential")
    test_model1 = [1 if abs(potential[x[0]]-potential[x[1]]) > 0.05 else 0 for x in node_index]

    electrons = get_node_model_values(device=device, region=region, name="Electrons")
    test_model2 = [1 if abs(math.log10(electrons[x[0]])-math.log10(electrons[x[1]])) > 1 else 0 for x in node_index]

    test_model = max_merge_lists((test_model1, test_model2))

    return test_model

def max_merge_lists(list_of_lists):
    '''
    gets the max of list of lists
    '''
    test_model = [max(x) for x in zip(*list_of_lists)]
    return test_model

def refine_silicon_region(fh, device, region, mincl, maxcl):
    '''
    refinement for silicon regions
    mincl : minimum characteristic length
    maxcl : maximum characteristic length
    '''

    test_model = get_silicon_model_values(device=device, region=region)

    refine_common(fh=fh, device=device, region=region, model_values=test_model, mincl=mincl, maxcl=maxcl)



