#!/usr/bin/env python
import json

'''
Created on Jun 4, 2015

@author: dlemmer
'''

ASSAY_TYPES = ("presence/absence", "SNP", "gene variant", "ROI", "mixed")

TARGET_FUNCTIONS = ("species ID", "strain ID", "resistance type", "virulence factor")

class Assay(object):
    '''
    classdocs
    '''

    def __init__(self, name, assay_type, target=None, percid=None):
        '''
        Constructor
        '''
        self.name = name
        self.assay_type = assay_type
        self.target = target
        self.percid = percid

    def __str__(self):
        return "Assay: %s (%s) \n{\n%s}\n" % (self.name, self.assay_type, self.target)

    @property
    def assay_type(self):
        return self._assay_type

    @assay_type.setter
    def assay_type(self, value):
        if value and value not in ASSAY_TYPES:
            raise TypeError('Not a valid Assay Type')
        self._assay_type = value

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, value):
        if value and not isinstance(value, Target):
            raise TypeError('Not a valid Target')
        self._target = value

    @property
    def percid(self):
        return self._percid

    @percid.setter
    def percid(self, value):
        if value and type(value) != float:
            raise TypeError('Not a valid percent identity')
        self._percid = value

class Target(object):
    '''
    classdocs
    '''

    def __init__(self, function, gene_name=None, start_position=None, end_position=None, reverse_comp=False, amplicon=None):
        '''
        Constructor
        '''
        self.function = function
        self.gene_name = gene_name
        self.start_position = start_position
        self.end_position = end_position
        self.reverse_comp = reverse_comp
        self.amplicons = []
        if amplicon:
            if isinstance(amplicon, list):
                for amp in amplicon:
                    self.add_amplicon(amp)
            else:
                self.add_amplicon(amplicon)

    def __str__(self):
        output = "Target: %s" % self.function
        if self.gene_name:
            output += " (%s" % self.gene_name
            if self.start_position:
                direction = "<-" if self.reverse_comp else "->"
                output += " (:%d-%d %s" % (self.start_position, self.end_position, direction)
            output += ")"
        for amplicon in self.amplicons:
            output += "\n\t%s" % amplicon
        return output

    def as_dict(self):
        return_dict = {'function':self._function}
        if self.gene_name:
            return_dict['gene_name'] = self.gene_name
        if self.start_position:
            return_dict['start_position'] = self.start_position
        if self.end_position:
            return_dict['end_position'] = self.end_position
        if self.reverse_comp:
            return_dict['reverse_comp'] = True
        return return_dict

    @property
    def function(self):
        return self._function

    @function.setter
    def function(self, value):
        if value and value not in TARGET_FUNCTIONS:
            raise TypeError('Not a valid Target Function')
        self._function = value

    def add_amplicon(self, amplicon):
        if not isinstance(amplicon, Amplicon):
            raise TypeError('Not a valid Amplicon')
        self.amplicons.append(amplicon)

class Amplicon(object):
    '''
    classdocs
    '''

    def __init__(self, sequence, variant_name=None, significance=None, snp=None, regionofinterest=None, percid=None):
        '''
        Constructor
        '''
        self.sequence = sequence
        self.variant_name = variant_name
        self.significance = significance
        self.percid = percid
        self.SNPs = []
        self.ROIs = []
        if snp:
            if isinstance(snp, list):
                for each_snp in snp:
                    self.add_SNP(each_snp)
            else:
                self.add_SNP(snp)
        if regionofinterest:
            if isinstance(regionofinterest, list):
                for roi in regionofinterest:
                    self.add_ROI(roi)
            else:
                self.add_ROI(regionofinterest)

    def __str__(self):
        output = "Amplicon: "
        if self.variant_name:
            output += "%s = %s\n" % (self.variant_name, self.significance)
        else:
            output += "{"
            output += self.sequence+" "
            for snp in self.SNPs:
                output += "%s, " % snp
            for roi in self.ROIs:
                output += "%s, " % roi
            if self.significance:
                output += "= %s" % self.significance
            output += "}\n"
        return output

    def as_dict(self):
        return_dict = {'sequence':self.sequence}
        if self.variant_name:
            return_dict['variant_name'] = self.variant_name
        if self._significance:
            return_dict['Significance'] = self._significance
        if self.percid:
            return_dict['percid'] = self.percid
        return return_dict

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if value:
            value = value.replace('\n', '')
            value = value.replace(' ', '')
            value = value.upper()
        self._sequence = value

    @property
    def significance(self):
        return self._significance

    @significance.setter
    def significance(self, value):
        if value and not isinstance(value, Significance):
            raise TypeError('Not a valid Significance')
        self._significance = value

    @property
    def percid(self):
        return self._percid

    @percid.setter
    def percid(self, value):
        if value and type(value) != float:
            raise TypeError('Not a valid percent identity')
        self._percid = value

    def add_SNP(self, snp):
        if not isinstance(snp, SNP):
            raise TypeError('Not a valid SNP')
        self.SNPs.append(snp)

    def add_ROI(self, roi):
        if not isinstance(roi, RegionOfInterest):
            raise TypeError('Not a valid RegionOfInterest')
        self.ROIs.append(roi)

class SNP(object):
    '''
    classdocs
    '''

    def __init__(self, position, reference, variant=None, name=None, significance=None):
        '''
        Constructor
        '''
        self.position = int(position)
        self.reference = reference
        self.variant = variant
        self.name = name
        self.significance = significance

    def __str__(self):
        return "SNP: %d%s->%s = %s" % (self.position, self.reference, self.variant, self.significance)

    def as_dict(self):
        return_dict = {'position':self.position, 'reference':self.reference}
        if self.variant:
            return_dict['variant'] = self.variant
        if self.name:
            return_dict['name'] = self.name
        if self._significance:
            return_dict['Significance'] = self._significance
        return return_dict

    @property
    def significance(self):
        return self._significance

    @significance.setter
    def significance(self, value):
        if value and not isinstance(value, Significance):
            raise TypeError('Not a valid Significance')
        self._significance = value

class RegionOfInterest(object):
    '''
    classdocs
    '''

    def __init__(self, position_range, aa_sequence=None, nt_sequence=None, mutations=None, name = None, significance=None):
        '''
        Constructor
        '''
        self.position_range = position_range
        self.aa_sequence = aa_sequence
        self.nt_sequence = nt_sequence
        self.mutations = mutations
        self.name = name
        self.significance = significance

    def as_dict(self):
        return_dict = {'position_range':self.position_range}
        if self.aa_sequence:
            return_dict['aa_sequence'] = self.aa_sequence
        if self.nt_sequence:
            return_dict['nt_sequence'] = self.nt_sequence
        if self.mutations:
            return_dict['mutations'] = ','.join(self.mutations)
        if self.name:
            return_dict['name'] = self.name
        if self._significance:
            return_dict['Significance'] = self._significance
        return return_dict

    @property
    def mutations(self):
        return self._mutations

    @mutations.setter
    def mutations(self, value):
        import re
        if not value or re.match('any', value, re.IGNORECASE):
            value = []
        else:
            value = value.split(',')
        self._mutations = value

    @property
    def significance(self):
        return self._significance

    @significance.setter
    def significance(self, value):
        if value and not isinstance(value, Significance):
            raise TypeError('Not a valid Significance')
        self._significance = value

class Significance(object):
    '''
    classdocs
    '''

    def __init__(self, message, resistance=None):
        '''
        Constructor
        '''
        self.message = message
        self.resistance = resistance

    def __str__(self):
        return "Significance: %s" % self.message

    def as_dict(self):
        return_dict = {'message':self.message}
        if self.resistance:
            return_dict['resistance'] = self.resistance
        return return_dict

# class AND(object):
#     '''
#     classdocs
#     '''
#
#     def __init__(self, significance, name=None, target=None, snp=None, regionofinterest=None):
#         '''
#         Constructor
#         '''
#         self.significance = significance
#         self.name = name
#         self.operands = []
#         if target:
#             if isinstance(target, list):
#                 for each_target in target:
#                     self.add_operand(each_target)
#             else:
#                 self.add_operand(target)
#         if snp:
#             if isinstance(snp, list):
#                 for each_snp in snp:
#                     self.add_operand(each_snp)
#             else:
#                 self.add_operand(snp)
#         if regionofinterest:
#             if isinstance(regionofinterest, list):
#                 for roi in regionofinterest:
#                     self.add_operand(roi)
#             else:
#                 self.add_operand(regionofinterest)
#
#     def __str__(self):
#         return "AND: {" + ", ".join(map(str, self.operands)) + "} = %s" % self.significance
#
#     @property
#     def significance(self):
#         return self._significance
#
#     @significance.setter
#     def significance(self, value):
#         if value and not isinstance(value, Significance):
#             raise TypeError('Not a valid Significance')
#         self._significance = value
#
#     def add_operand(self, value):
#         if not (isinstance(value, Target) or isinstance(value, SNP) or isinstance(value, RegionOfInterest)):
#             raise TypeError('Not a valid operand')
#         self.operands.append(value)


class Operation(object):
    '''
    A representation of an operation, such as AND, OR, NOT, etc
    'type' represents the operand ie AND
    'children' represent a list of either a nested operation or a item to evaluate
    Example: a AND b would be a operation of type 'AND' with children a, b as ITEMs
    '''

    def __init__(self, operation_type, children, message=""):
        knownOperations = ["AND", "OR", "NOT"]
        '''
        Constructor
        Required:
            type=<String of any of these options %s>,
            children=[<ITEM or Operation>]
        ''' % str(knownOperations)
        self.message = message
        if operation_type:
            operation_type = str(operation_type).upper()
            pass
        if not isinstance(operation_type, str):
            raise AttributeError("Operation type must be a string")
            pass
        if operation_type not in knownOperations:
            raise AttributeError("Operation type must be of a known operation: %s" % str(knownOperations))
            pass
        self.operation_type = operation_type
        if not isinstance(children, list):
            raise AttributeError("Operation children must be a list")
            pass
        self.children = []
        for child in children:
            if "item" in child:
                if isinstance(child["item"], list):
                    childItemList = child["item"]
                    for childItem in childItemList:
                        if not isinstance(childItem, ITEM):
                            AttributeError("Operation item child was not parsed as ITEM type")
                            pass
                        self.children.append(childItem)
                        pass
                    pass
                else:
                    if not isinstance(child["item"], ITEM):
                        AttributeError("Operation item child was not parsed as ITEM type")
                        pass
                    self.children.append(child["item"])

                pass
            if "operation" in child:
                if isinstance(child["operation"], list):
                    childOperList = child["operation"]
                    for childOper in childOperList:
                        if not isinstance(childOper, Operation):
                            AttributeError("Operation child was not a proper operation")
                            pass
                        self.children.append(childOper)
                        pass
                else:
                    if not isinstance(child["operation"], Operation):
                        AttributeError("Operation child was not a proper operation")
                        pass
                    self.children.append(child["operation"])
                pass
            if isinstance(child, Operation) or isinstance(child, ITEM):
                self.children.append(child) # This should always be a list of ITEMs or Operations
                pass
            pass

        if operation_type == "NOT":
            if len(self.children) > 1:
                raise AttributeError("Operation type of 'NOT' can only have 1 child")
                pass
            pass
        else:
            if len(self.children) < 2:
                raise AttributeError("Operation type needs at least 2 children")
                pass
            pass



    def __str__(self):
        childStr = ""
        for child in self.children:
            if len(childStr) < 1:
                childStr = str(child)
                pass
            else:
                childStr = childStr+" , "+str(child)
            pass
        return "Operation: ['operation_type':%s, children:[%s]]" % (self.operation_type, childStr)


class ITEM(object):
    '''
    A representation of a known data type within asap as well as a key and value to find it by
        ie key="name" value="Ba_specific-6_CP012725.1_3629" means to look for "name:Ba_specific-6_CP012725.1_3629"
        key can also define how the item is evaluted ie key="depth_greater_than" will be treated as true
        if the depth value in the data is greater than the value providied for the Item
    Used in creating operation dictionaries
    '''

    def __init__(self, item_type, identity_key, identity_value, evaluation=None, value=None):
        '''
        Constructor
        '''
        # dataType, and identity variables are used to find the entry
        self.item_type = str(item_type).lower()
        self.identity_key = str(identity_key).lower()
        self.identity_value = str(identity_value).lower()
        # evaluation and value are used to add a new check beyond ensuring existiance of data
        self.evaluation = str(evaluation).lower()
        self.value = str(value).lower()

    def __str__(self):
        return "ITEM:{item_type=%s, identity_key=%s, identity_value=%s, evaluation=%s, value=%s}" % (str(self.item_type), str(self.identity_key), str(self.identity_value), str(self.evaluation), str(self.value))


def _json_decode(json_dict):
    #print("json_dict = %s" % json_dict)
    json_dict = dict((k.lower() if k != 'AND' else 'AND', v) for k,v in json_dict.items())
    if "message" in json_dict:
        return Significance(**json_dict)
    elif "position_range" in json_dict:
        return RegionOfInterest(**json_dict)
    elif "position" in json_dict:
        return SNP(**json_dict)
    elif "sequence" in json_dict:
        return Amplicon(**json_dict)
    elif "function" in json_dict:
        return Target(**json_dict)
    elif "assay_type" in json_dict:
        return Assay(**json_dict)
    elif "operation_type" in json_dict:
        # print("OPERATION: ", json_dict)
        return Operation(**json_dict)
    elif "item_type" in json_dict:
        # print("ITEM: ", json_dict)
        return ITEM(**json_dict)
    else:
        return json_dict

# This needs to be beefed up a lot, but fine for converting fasta to json
def _json_encode(obj):
    if isinstance(obj, Assay):
        assay_dict = {"name":obj.name, "assay_type":obj.assay_type}
        if obj.target:
            assay_dict["Target"] = obj.target
        if obj.percid:
            assay_dict["percid"] = obj.percid
        return assay_dict
    if isinstance(obj, Target):
        target_dict = obj.as_dict()
        if obj.amplicons:
            target_dict["Amplicon"] = obj.amplicons if len(obj.amplicons) > 1 else obj.amplicons[0]
        return target_dict
    if isinstance(obj, Amplicon):
        amplicon_dict = obj.as_dict()
        if obj.SNPs:
            amplicon_dict["SNP"] = obj.SNPs if len(obj.SNPs) > 1 else obj.SNPs[0]
        if obj.ROIs:
            amplicon_dict["RegionOfInterest"] = obj.ROIs if len(obj.ROIs) > 1 else obj.ROIs[0]
        if obj.percid:
            amplicon_dict["percid"] = obj.percid
        return amplicon_dict
    if isinstance(obj, SNP):
        snp_dict = obj.as_dict()
        return snp_dict
    if isinstance(obj, RegionOfInterest):
        roi_dict = obj.as_dict()
        return roi_dict
    if isinstance(obj, Significance):
        sig_dict = obj.as_dict()
        return sig_dict
    else:
        return json.JSONEncoder.default(obj)

def parseJSON(file_path):
    with open(file_path, 'r') as file_object:
        return json.load(file_object, object_hook=_json_decode)['assay']

def parseOperation(file_path):
    with open(file_path, 'r') as file_object:
        return json.load(file_object, object_hook=_json_decode)['operation']

def writeJSON(assay_data, filename):
    with open(filename, 'w') as json_fh:
        json.dump(assay_data, json_fh, indent=2, default=_json_encode)

def generateReference(assay_list):
    from skbio import DNA
    for assay in assay_list:
        name = assay.name
        for amplicon in assay.target.amplicons:
            name = assay.name + "_%s" % amplicon.variant_name if amplicon.variant_name else name
            seq = DNA(amplicon.sequence, {'id':name})
            yield seq


def main():
    import argparse
    import json
    import skbio
    import sys
    parser = argparse.ArgumentParser(description='backdoor command to generate a reference fasta from an ASAP JSON file', epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('json', type=argparse.FileType('r'), help='JSON file to print as a fasta')
    args = parser.parse_args()

    data = json.load(args.json, object_hook=_json_decode)
    skbio.io.registry.write(generateReference(data['assay']), "fasta", sys.stdout)

if __name__ == "__main__":
    main()
