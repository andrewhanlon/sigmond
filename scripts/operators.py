import xml.etree.ElementTree as ET
import xmldict

class Operator:

  def __init__(self, op_xml):

    self.op_xml = op_xml
    self.op_xml_strip = ET.tostring(op_xml).decode().replace('\n', '').replace(' ', '')
    self.op_string = self.__xml_to_string(op_xml)

  def __xml_to_string(self, op_xml):

    num_hadrons = int(op_xml.find("NumberOfHadrons").text)

    op_string = ''

    if (num_hadrons == 1):
      hadron = xmldict.xml_to_dict(ET.tostring(op_xml.find("Hadron")).decode())['Hadron']
      irrep_row = op_xml.find("LGIrrepRow").text
      mom = hadron['Momentum'].split()

      spatial = ''
      if (hadron['Flavor'] == "glueball"):
        spatial = hadron['SpatialType']
      else:
        spatial = hadron['SpatialType'] + "_" + hadron['SpatialIdNum']

      op_string = hadron['Flavor'] + " P=(" + mom[0] + "," + mom[1] + "," + mom[2] + ") " + hadron['LGIrrep'] + \
                  "_" + irrep_row + " " + spatial
      

    elif (num_hadrons == 2):

      total = xmldict.xml_to_dict(ET.tostring(op_xml.find("Total")).decode())['Total']
      hadron1 = xmldict.xml_to_dict(ET.tostring(op_xml.find("Hadron1")).decode())['Hadron1']
      hadron2 = xmldict.xml_to_dict(ET.tostring(op_xml.find("Hadron2")).decode())['Hadron2']

      had1_mom = hadron1['Momentum'].split()
      had2_mom = hadron2['Momentum'].split()

      lg_cg_id = ""
      if 'LGCGId' in total:
        lg_cg_id = " CG_" + total['LGCGId']

      op_string = "iso" + total['Isospin'] + "_" + hadron1['Flavor'] + "_" + hadron2['Flavor'] + " " + \
                  total['LGIrrep'] + "_" + total['LGIrrepRow'] + lg_cg_id + \
                  " [P=(" + had1_mom[0] + "," + had1_mom[1] + "," + had1_mom[2] + ") " + hadron1['LGIrrep'] + " " + \
                  hadron1['SpatialType'] + "_" + hadron1['SpatialIdNum'] + \
                  "] [P=(" + had2_mom[0] + "," + had2_mom[1] + "," + had2_mom[2] + ") " + hadron2['LGIrrep'] + " " + \
                  hadron2['SpatialType'] + "_" + hadron2['SpatialIdNum'] + "]"


    else:

      print("WARNING: Operator has neither 1 nor 2 hadrons")

    return op_string


  def __repr__(self):

    return self.op_string

  def __eq__(self, other):

    if isinstance(other, Operator):
      return (self.op_string == other.op_string)
    else:
      return False

  def __ne__(self, other):
    
    return (not self.__eq__(other))

  def __hash__(self):

    return hash(self.__repr__())
