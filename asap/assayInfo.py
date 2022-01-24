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

    def __init__(self, name, assay_type, target=None, AND=None):
        '''
        Constructor
        '''
        self.name = name
        self.assay_type = assay_type
        self.target = target
        self.AND = AND
    
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
    def AND(self):
        return self._AND
    
    @AND.setter        
    def AND(self, value):
        if value and not isinstance(value, AND):
            raise TypeError('Not a valid AND operation')
        self._AND = value
        
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

    def __init__(self, sequence, variant_name=None, significance=None, snp=None, regionofinterest=None, AND=None):
        '''
        Constructor
        '''
        self.sequence = sequence
        self.variant_name = variant_name
        self.significance = significance
        self.SNPs = []
        self.ROIs = []
        self.AND = AND
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
            if self.AND:
                output += "%s" % self.AND
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
        
    def add_SNP(self, snp):
        if not isinstance(snp, SNP):
            raise TypeError('Not a valid SNP')
        self.SNPs.append(snp)
        
    def add_ROI(self, roi):
        if not isinstance(roi, RegionOfInterest):
            raise TypeError('Not a valid RegionOfInterest')
        self.ROIs.append(roi)
        
    @property
    def AND(self):
        return self._AND
    
    @AND.setter        
    def AND(self, value):
        if value and not isinstance(value, AND):
            raise TypeError('Not a valid AND operation')
        self._AND = value
        
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
    
class AND(object):
    '''
    classdocs
    '''

    def __init__(self, significance, name=None, target=None, snp=None, regionofinterest=None):
        '''
        Constructor
        '''
        self.significance = significance
        self.name = name
        self.operands = []
        if target:
            if isinstance(target, list):
                for each_target in target:
                    self.add_operand(each_target)
            else:
                self.add_operand(target)
        if snp:
            if isinstance(snp, list):
                for each_snp in snp:
                    self.add_operand(each_snp)
            else:
                self.add_operand(snp)
        if regionofinterest:
            if isinstance(regionofinterest, list):
                for roi in regionofinterest:
                    self.add_operand(roi)
            else:
                self.add_operand(regionofinterest)
    
    def __str__(self):
        return "AND: {" + ", ".join(map(str, self.operands)) + "} = %s" % self.significance          
        
    @property
    def significance(self):
        return self._significance
    
    @significance.setter
    def significance(self, value):
        if value and not isinstance(value, Significance):
            raise TypeError('Not a valid Significance')
        self._significance = value

    def add_operand(self, value):
        if not (isinstance(value, Target) or isinstance(value, SNP) or isinstance(value, RegionOfInterest)):
            raise TypeError('Not a valid operand')
        self.operands.append(value)
        

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
    elif "snp" in json_dict or "regionofinterest" in json_dict or "target" in json_dict:
        return AND(**json_dict)
    else:
        return json_dict

# This needs to be beefed up a lot, but fine for converting fasta to json    
def _json_encode(obj):
    if isinstance(obj, Assay):
        assay_dict = {"name":obj.name, "assay_type":obj.assay_type}
        if obj.target:
            assay_dict["Target"] = obj.target
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
    
def parseJSON(file_object):
    return json.load(file_object, object_hook=_json_decode)['assay']
    
def writeJSON(assay_data, filename):
    with open(filename, 'w') as json_fh:
        json.dump(assay_data, json_fh, indent=2, default=_json_encode)
    
def generateReference(assay_list):
    from skbio import DNA
    for assay in assay_list:
        name = assay.name
        if assay.AND:
            for operand in assay.AND:
                if isinstance(operand, Target):
                    name = name + "_%s" % operand.gene_name if operand.gene_name else name
                    for amplicon in operand.amplicons:
                        name = name + "_%s" % amplicon.variant_name if amplicon.variant_name else name
                        seq = DNA(amplicon.sequence, id=name)
                        yield seq
        else:
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
