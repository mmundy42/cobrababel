"""
@author: Jay
"""
import re
import warnings


def translate(model, reaction_xref_file_name, metabolite_xref_file_name, from_namespace, to_namespace):
    
    """ Translate IDs in a model from one namespace to another namespace

    Parameters
    ----------
    model : cobra.core.Model
        Model to translate
    reaction_xref_file_name : str
        Path to cross reference file with ID mapping for reactions
    metabolite_xref_file_name : str
        Path to cross reference file with ID mapping for metabolites
    from_namespace : str
        Namespace of IDs in input model
    to_namespace : str
        Namespace of IDs in output model

    Returns
    -------
    cobra.core.Model
        New model with translated IDs
    """   
    
    ''' reading input cobra model ''' 
    #model="C:/Users/M179100/GitHub/cobrababel/data/Btheta.xml" 
    output_model = model
    '''
    ######################################################    
            converting name space for metabolites     
    ######################################################
    '''
    metabolite_conversion_lines = open(metabolite_xref_file_name).readlines()
    metabolite_conversion_header = metabolite_conversion_lines[0].strip().split("\t")
    
    ''' raise exception for invalid name space files '''
    if len(metabolite_conversion_header)!= 2:
        raise ValueError('metabolite name space reference file is invalid')    
    metabolite_conversion_lines = metabolite_conversion_lines[1:]
    
    ''' creating name space dictionary of metabolites from cross reference file '''
    metabolite_seed2vmh_dict = {}
    if from_namespace == metabolite_conversion_header[0] and to_namespace == metabolite_conversion_header[1]:
        for m in metabolite_conversion_lines:
            s,v = m.split("\t")
            s,v = s.strip("\n"),v.strip("\n")
            #cs,cv=re.sub("[EX_ (e)]","",s).strip(),re.sub("[EX_ (e)]","",v).strip()  
            metabolite_seed2vmh_dict[s] = v
    elif from_namespace == metabolite_conversion_header[1] and to_namespace == metabolite_conversion_header[0]:
         for m in metabolite_conversion_lines:
            s,v = m.split("\t")
            s,v = s.strip("\n"),v.strip("\n")
            #cs,cv=re.sub("[EX_ (e)]","",s).strip(),re.sub("[EX_ (e)]","",v).strip()  
            metabolite_seed2vmh_dict[v] = s
    else :
            raise ValueError("from and to name spaces are invalid for metabolites, conversion not possible")
    
    ''' replacing the name space ids for metabolites using  dictionary '''
    out_metabolites = output_model.metabolites
    suffix = re.compile(r'_([ce])$')        
    metabolites_not_found_in_Xref = []
    metabolites_found_in_Xref = []
    metabolites_already_converted = []
    for o_m in out_metabolites:
        try:
            #o_m=out_metabolites[0]
            s = o_m.id 
            me,su,ex = re.split(suffix,s)
            if me in metabolite_seed2vmh_dict.keys():
                v = metabolite_seed2vmh_dict[me].strip()
                o_m.id = v + "_" + su
                metabolites_found_in_Xref.append(str(me + " : " + o_m.id))
                #print s,v
            else:
                if me not in metabolite_seed2vmh_dict.values():
                    metabolites_not_found_in_Xref.append(me)
                else : 
                    metabolites_already_converted.append(me)
        except:
            pass
    
    ''' A warning for non converted metabolites count '''
    warn_text_m = "Could not convert name space for "+str(len(metabolites_not_found_in_Xref))+" metabolites"
    warnings.warn(warn_text_m)
    
    
    '''
    ###########################################################    
                converting namespace for reactions
    ###########################################################
    '''
    reaction_coversion_lines = open(reaction_xref_file_name).readlines()
    reaction_conversion_header = reaction_coversion_lines[0].strip().split("\t")
    
    '''raise exception  for invalid name space files '''
    if len(reaction_conversion_header)!= 2:
        raise ValueError('reaction name space reference file is invalid')
    reaction_coversion_lines = reaction_coversion_lines[1:]
    
    ''' creating name space dictionary of reactions from cross reference file '''
    reactions_seed2vmh_dict = {}
    if from_namespace == reaction_conversion_header[0] and to_namespace == reaction_conversion_header[1]:
        for r in reaction_coversion_lines:
            s,v = r.split("\t")
            cs,cv = s.strip(),v.strip()
            #cs,cv=s.replace("_e","").strip(),v.replace("_e","").strip()
            #cs,cv=re.sub("[(e)]","_e",s).strip(),re.sub("[(e)]","_e",v).strip()  
            reactions_seed2vmh_dict[cs] = cv
    elif from_namespace == reaction_conversion_header[1] and to_namespace == reaction_conversion_header[0]:
         for r in reaction_coversion_lines:
            s,v = r.split("\t")
            cs,cv = s.strip(),v.strip()
            #cs,cv=s.replace("_e","").strip(),v.replace("_e","").strip()  
            reactions_seed2vmh_dict[cv] = cs
    else :
         raise ValueError("from and to name spaces are invalid for reactions, conversion not possible")
    
    
    ''' replacing the name space ids for reactions using  dictionary '''
    out_reactions = output_model.reactions
    reactions_not_found_in_Xref = []
    reactions_found_in_Xref = []
    reactions_already_converted = []
    for o_r  in out_reactions:
        try:
            #o_r=out_reactions[0]
            s = o_r.id
            rea = s
            #rea,su,ex=re.split(suffix,s)
            if rea in reactions_seed2vmh_dict.keys():
                v = reactions_seed2vmh_dict[rea].strip()
                o_r.id = v #+"_"+su
                reactions_found_in_Xref.append(str(rea +"  :  " +o_r.id))
                #print rea,o_r.id
            else:
                if rea not in reactions_seed2vmh_dict.values():
                    reactions_not_found_in_Xref.append(rea)
                else:
                    reactions_already_converted.append(rea)
        except:
            pass 

    warn_text_m = "Could not convert name space for "+str(len(reactions_not_found_in_Xref))+" reactions"
    warnings.warn(warn_text_m)
    return output_model
