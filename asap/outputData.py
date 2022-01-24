'''
Created on Jun 17, 2015

@author: dlemmer
'''

import logging
from xml.etree import ElementTree


def _write_parameters(node, data):
    for k, v in data.items():
        subnode = ElementTree.SubElement(node, k)
        subnode.text = v
    return node

def _parse_parameters(node):
    parms = {}
    for element in node.iter():
        parms[element.tag] = element.text
    return parms



def writeOutput(parms, xml_file):
    from xml.dom import minidom

    root = ElementTree.Element("matrix_data")
    parm_node = ElementTree.SubElement(root, "parameters")
    _write_parameters(parm_node, parms)
    files_node = ElementTree.SubElement(root, "files")
    dom = minidom.parseString(ElementTree.tostring(root, 'utf-8'))
    output = open(xml_file, 'w')
    output.write(dom.toprettyxml(indent="    "))
    output.close()
    return xml_file

def createXMLNode(name, attributes):
    return ElementTree.Element(name, attributes)


def parse(xml_file):
    xmltree = ElementTree.parse(xml_file)
    root = xmltree.getroot()
    parms = _parse_parameters(root.find("parameters"))
    return parms

if __name__ == "__main__":
    pass
