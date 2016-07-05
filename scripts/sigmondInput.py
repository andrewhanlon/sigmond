import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import os
import subprocess
import operators

class SigmondInput:

  def __init__(self, projName, direc, echo, bootstrap, laph_query):

    self.ops = set()
    self.has_vevs = False
    self.herm_corr = False
    self.direc = direc
    self.channel_name = self.__channel_name()
    self.proj_name = projName
    self.full_name = self.channel_name + "_" + self.proj_name

    root = ET.Element("SigMonD")
    init_tag = ET.SubElement(root, "Initialize")
    proj_name_tag = ET.SubElement(init_tag, "ProjectName")
    proj_name_tag.text = self.full_name
    log_tag = ET.SubElement(init_tag, "LogFile")
    log_tag.text = self.full_name + ".log"

    if echo:
      ET.SubElement(init_tag, "EchoXML")

    ET.SubElement(init_tag, "MCObservables")

    self.bootstrap = bootstrap
    if self.bootstrap:
      boot_tag = ET.SubElement(init_tag, "Bootstrapper")
      num_samplings_tag = ET.SubElement(boot_tag, "NumberResamplings")
      num_samplings_tag.text = "2048"
      seed_tag = ET.SubElement(boot_tag, "Seed")
      seed_tag.text = "6754"
      skip_tag = ET.SubElement(boot_tag, "BootSkip")
      skip_tag.text = "37"
      ET.SubElement(boot_tag, "Precompute")

    ET.SubElement(root, "TaskSequence")
    self.xml_root = root
    
    self.add_data_files(laph_query)

  def __add__(self, sig2):
    
    sig1 = copy.deepcopy(self)

    for fileListInfo in sig2.xml_root.find("CorrelatorData"):
      sig1.xml_root.find("CorrelatorData").append(fileListInfo)

    for fileListInfo in sig2.xml_root.find("VEVData"):
      sig1.xml_root.find("VEVData").append(fileListInfo)

    for task in sig2.xml_root.find("TaskSequence"):
      sig1.xml_root.find("TaskSequence").append(task)

    return sig1


  def add_data_files(self, laph_query):

    corr_files = []
    vev_files = []

    for f in os.listdir(self.direc):
      
      full_file = os.path.join(self.direc, f)

      if not os.path.isfile(full_file):
        continue

      if self.__is_corr(full_file):
        corr_files.append(full_file)
      elif self.__is_vev(full_file):
        vev_files.append(full_file)

        
    if len(corr_files):
      self.__add_data_files(corr_files, True)

    if len(vev_files):
      self.has_vevs = True
      self.__add_data_files(vev_files, False)


    self.__set_operators(corr_files, vev_files, laph_query)

  def print_ops(self):
    
    print(self.channel_name, ":\n")
    for op in self.ops:
      print(op.op_string)

    print("\n\n")

  def print_ops(self, filename):

    f = open(filename, 'w')
    
    for op in self.ops:
      f.write(op.op_string + '\n')

    f.close()

  def keep_ops(self, filename):

    f = open(filename, 'r')

    keep_ops = set()

    for line in f:
      keep_ops.add(line.rstrip())

    if '' in keep_ops:
      keep_ops.remove('')

    temp_ops = self.ops.copy()

    for op in temp_ops:
      if op.op_string not in keep_ops:
        self.ops.remove(op)

    f.close()
      
  def write(self):
    
    xml_file = self.channel_name + "_" + self.proj_name + ".xml"

    xmlstr = minidom.parseString(ET.tostring(self.xml_root)).toprettyxml(indent="  ")
    with open(xml_file, "w") as f:
      f.write(xmlstr)

    print("Saved input xml to: " + xml_file)

  def execute(self, sigmond):
    
    xml_file = self.channel_name + "_" + self.proj_name + ".xml"
    print("Executing Sigmond...")
    subprocess.call(['sigmond', xml_file])

  def read(self):

    print("Reading results...\n")
      

  def do_checks(self):

    task_seq_tag = self.xml_root.find("TaskSequence")
    task_tag = ET.SubElement(task_seq_tag, "Task")
    action_tag = ET.SubElement(task_tag, "Action")
    action_tag.text = "DoChecks"
    type_tag = ET.SubElement(task_tag, "Type")
    type_tag.text = "TemporalCorrelatorMatrix"
    corr_tag = ET.SubElement(task_tag, "CorrelatorMatrixInfo")

    if self.herm_corr:
      ET.SubElement(corr_tag, "HermitianMatrix")

      if self.has_vevs:
        ET.SubElement(corr_tag, "SubtractVEV")

    for op in self.ops:

      op_tag = ET.SubElement(corr_tag, "OperatorString")
      op_tag.text = op.op_string

    min_time_sep_tag = ET.SubElement(task_tag, "MinTimeSep")
    min_time_sep_tag.text = "3"
    max_time_sep_tag = ET.SubElement(task_tag, "MaxTimeSep")
    max_time_sep_tag.text = "25"
    verbose_tag = ET.SubElement(task_tag, "Verbose")
    outlier_tag = ET.SubElement(task_tag, "OutlierScale")
    outlier_tag.text = "15"
    

    task_tag = ET.SubElement(task_seq_tag, "Task")
    action_tag = ET.SubElement(task_tag, "Action")
    action_tag.text = "DoChecks"
    type_tag = ET.SubElement(task_tag, "Type")
    type_tag.text = "TemporalCorrelatorMatrixIsHermitian"
    corr_tag = ET.SubElement(task_tag, "CorrelatorMatrixInfo")

    for op in self.ops:
      
      op_tag = ET.SubElement(corr_tag, "OperatorString")
      op_tag.text = op.op_string

    min_time_sep_tag = ET.SubElement(task_tag, "MinTimeSep")
    min_time_sep_tag.text = "3"
    max_time_sep_tag = ET.SubElement(task_tag, "MaxTimeSep")
    max_time_sep_tag.text = "25"
    verbose_tag = ET.SubElement(task_tag, "Verbose")

  
  def set_herm(self):
    
    self.herm_corr = True

  def do_fits(self, plot):

    idn = 0
    plot_file_num = 0

    task_seq_tag = self.xml_root.find("TaskSequence")

    for op in self.ops:

      task_tag = ET.SubElement(task_seq_tag, "Task")
      action_tag = ET.SubElement(task_tag, "Action")
      action_tag.text = "DoFit"

      type_tag = ET.SubElement(task_tag, "Type")
      type_tag.text = "TemporalCorrelator"

      min_info_tag = ET.SubElement(task_tag, "MinimizerInfo")
      min_method_tag = ET.SubElement(min_info_tag, "Method")
      min_method_tag.text = "Minuit2"
      min_param_tag = ET.SubElement(min_info_tag, "ParameterRelTol")
      min_param_tag.text = "2e-6"
      min_chi_tag = ET.SubElement(min_info_tag, "ChiSquareRelTol")
      min_chi_tag.text = "2e-4"
      min_iter_tag = ET.SubElement(min_info_tag, "MaximumIterations")
      min_iter_tag.text = "1000"
      min_verbose_tag = ET.SubElement(min_info_tag, "Verbosity")
      min_verbose_tag.text = "Low"

      sampling_tag = ET.SubElement(task_tag, "SamplingMode")
      if self.bootstrap:
        sampling_tag.text = "Boostrap"
      else:
        sampling_tag.text = "Jackknife"
        

      corr_fit_tag = ET.SubElement(task_tag, "TemporalCorrelatorFit")
      min_time_sep_tag = ET.SubElement(corr_fit_tag, "MinimumTimeSeparation")
      min_time_sep_tag.text = "10"
      max_time_sep_tag = ET.SubElement(corr_fit_tag, "MaximumTimeSeparation")
      max_time_sep_tag.text = "25"
      op_tag = ET.SubElement(corr_fit_tag, "OperatorString")
      op_tag.text = op.op_string
      noise_cutoff_tag = ET.SubElement(corr_fit_tag, "LargeTimeNoiseCutoff")
      noise_cutoff_tag.text = "1.0"
      model_tag = ET.SubElement(corr_fit_tag, "Model")
      model_type_tag = ET.SubElement(model_tag, "Type")
      model_type_tag.text = "TimeSymSingleExponentialPlusConstant"
      energy_tag = ET.SubElement(model_tag, "Energy")
      energy_name_tag = ET.SubElement(energy_tag, "Name")
      energy_name_tag.text = op.op_string.split()[0]
      energy_id_tag = ET.SubElement(energy_tag, "IDIndex")
      energy_id_tag.text = str(idn)
      idn += 1
      amp_tag = ET.SubElement(model_tag, "Amplitude")
      amp_name_tag = ET.SubElement(amp_tag, "Name")
      amp_name_tag.text = "A"
      amp_id_tag = ET.SubElement(amp_tag, "IDIndex")
      amp_id_tag.text = str(idn)
      idn += 1
      cons_tag = ET.SubElement(model_tag, "AddedConstant")
      cons_name_tag = ET.SubElement(cons_tag, "Name")
      cons_name_tag.text = "C0"
      cons_id_tag = ET.SubElement(cons_tag, "IDIndex")
      cons_id_tag.text = str(idn)
      idn += 1

      if plot:
        plot_tag = ET.SubElement(corr_fit_tag, "DoEffectiveEnergyPlot")
        plot_file_tag = ET.SubElement(plot_tag, "PlotFile")
        plot_file_tag.text = "Eff_energy" + str(plot_file_num) + "--" + op.op_string.split()[0] + ".agr"
        plot_file_num += 1
        energy_id_name_tag = ET.SubElement(plot_tag, "EffEnergyIdName")
        energy_id_name_tag.text = op.op_string.split()[0]
        time_sep_tag = ET.SubElement(plot_tag, "TimeStep")
        time_sep_tag.text = "1"
        sym_col_tag = ET.SubElement(plot_tag, "SymbolColor")
        sym_col_tag.text = "blue"
        sym_type_tag = ET.SubElement(plot_tag, "SybmolType")
        sym_type_tag.text = "circle"
        #plot_err_tag = ET.SubElement(plot_tag, "MaxErrorToPlot")
        #plot_err_tag.text = "

  def printXML(self, typePrint, real, bins, verbose):
    
    task_seq_tag = self.xml_root.find("TaskSequence")

    for op1 in self.ops:
        
      task_tag = ET.SubElement(task_seq_tag, "Task")
      action_tag = ET.SubElement(task_tag, "Action")
      action_tag.text = "PrintXML"
      type_tag = ET.SubElement(task_tag, "Type")
      type_tag.text = typePrint
      if typePrint.find("Histogram") >= 0:
        bins_tag = ET.SubElement(task_tag, "NumberOfBins")
        bins_tag.text = bins
      obs_tag = ET.SubElement(task_tag, "MCObservable")
      vev_tag = ET.SubElement(obs_tag, "VEV")
      src_op_tag = ET.SubElement(vev_tag, "OperatorString")
      src_op_tag.text = op1.op_string

      arg_tag = ET.SubElement(obs_tag, "Arg")
      if real:
        arg_tag.text = "Re"
      else:
        arg_tag.text = "Im"

      if verbose:
        ET.SubElement(task_tag, "Verbose")
      
      for op2 in self.ops:
        
        task_tag = ET.SubElement(task_seq_tag, "Task")
        action_tag = ET.SubElement(task_tag, "Action")
        action_tag.text = "PrintXML"
        type_tag = ET.SubElement(task_tag, "Type")
        type_tag.text = typePrint
        if typePrint.find("Histogram") >= 0:
          bins_tag = ET.SubElement(task_tag, "NumberOfBins")
          bins_tag.text = bins
        obs_tag = ET.SubElement(task_tag, "MCObservable")
        corr_tag = ET.SubElement(obs_tag, "Correlator")
        src_tag = ET.SubElement(corr_tag, "Source")
        src_op_tag = ET.SubElement(src_tag, "OperatorString")
        src_op_tag.text = op1.op_string
        snk_tag = ET.SubElement(corr_tag, "Sink")
        snk_op_tag = ET.SubElement(snk_tag, "OperatorString")
        snk_op_tag.text = op2.op_string
        time_tag = ET.SubElement(corr_tag, "TimeIndex")
        time_tag.text = "3"
        if self.herm_corr:
          ET.SubElement(corr_tag, "HermitianMatrix")
        if self.has_vevs:
          ET.SubElement(corr_tag, "SubtractVEV")


        arg_tag = ET.SubElement(obs_tag, "Arg")
        if real:
          arg_tag.text = "Re"
        else:
          arg_tag.text = "Im"

        if verbose:
          ET.SubElement(task_tag, "Verbose")

    

  def printEnergy(self, real):
    
    task_seq_tag = self.xml_root.find("TaskSequence")

    for op1 in self.ops:
      for op2 in self.ops:
        
        task_tag = ET.SubElement(task_seq_tag, "Task")
        action_tag = ET.SubElement(task_tag, "Action")
        action_tag.text = "PrintXML"
        type_tag = ET.SubElement(task_tag, "Type")
        type_tag.text = "EffectiveEnergy"
        eff_type_tag = ET.SubElement(task_tag, "EffEnergyType")
        eff_type_tag.text = "TimeForward" # @ADH - generalize this
        time_step_tag = ET.SubElement(task_tag, "TimeStep")
        time_step_tag.text = "3" # @ADH - generalize this
        # @ADH - implement this
        #id_tag = ET.SubElement(task_tag, "EffEnergyIdName")
        #id_tag.text = "TODO"
        corr_tag = ET.SubElement(task_tag, "Correlator")
        src_tag = ET.SubElement(corr_tag, "Source")
        src_op_tag = ET.SubElement(src_tag, "OperatorString")
        src_op_tag.text = op1.op_string
        snk_tag = ET.SubElement(corr_tag, "Sink")
        snk_op_tag = ET.SubElement(snk_tag, "OperatorString")
        snk_op_tag.text = op2.op_string

        arg_tag = ET.SubElement(task_tag, "Arg")
        if real:
          arg_tag.text = "Re"
        else:
          arg_tag.text = "Im"

        if self.herm_corr:
          ET.SubElement(task_tag, "HermitianMatrix")

        if self.has_vevs:
          ET.SubElement(task_tag, "SubtractVEV")

        sampling_tag = ET.SubElement(task_tag, "SamplingMode")
        if self.bootstrap:
          sampling_tag.text = "Bootstrap"
        else:
          sampling_tag.text = "Jackknife"
  
  
  def printCorr(self, real):

    task_seq_tag = self.xml_root.find("TaskSequence")

    for op1 in self.ops:
      for op2 in self.ops:
        
        task_tag = ET.SubElement(task_seq_tag, "Task")
        action_tag = ET.SubElement(task_tag, "Action")
        action_tag.text = "PrintXML"
        type_tag = ET.SubElement(task_tag, "Type")
        type_tag.text = "TemporalCorrelator"
        corr_tag = ET.SubElement(task_tag, "Correlator")
        src_tag = ET.SubElement(corr_tag, "Source")
        src_op_tag = ET.SubElement(src_tag, "OperatorString")
        src_op_tag.text = op1.op_string
        snk_tag = ET.SubElement(corr_tag, "Sink")
        snk_op_tag = ET.SubElement(snk_tag, "OperatorString")
        snk_op_tag.text = op2.op_string

        arg_tag = ET.SubElement(task_tag, "Arg")
        if real:
          arg_tag.text = "Re"
        else:
          arg_tag.text = "Im"

        if self.herm_corr:
          ET.SubElement(task_tag, "HermitianMatrix")

        if self.has_vevs:
          ET.SubElement(task_tag, "SubtractVEV")

        sampling_tag = ET.SubElement(task_tag, "SamplingMode")
        if self.bootstrap:
          sampling_tag.text = "Bootstrap"
        else:
          sampling_tag.text = "Jackknife"

 
  # 'Private' Methods
  
  def __channel_name(self):
    
    name = ''
    if "clover_s32" in self.direc:
      name += "32_"
    elif "clover_s24" in self.direc:
      name += "24_"
    else:
      print("WARNING: could not determine lattice ensemble")

    if "bosonic" in self.direc:
      name += "B_"
    elif "fermionic" in self.direc:
      name += "F_"
    else:
      print("WARNING: could not determine fermionic or bosonic")

    if "isosinglet" in self.direc:
      name += "I0_"
    elif "isodoublet" in self.direc:
      name += "2I1_"
    elif "isotriplet" in self.direc:
      name += "I1_"
    elif "isoquartet" in self.direc:
      name += "2I3_"
    elif "isoquintet" in self.direc:
      name += "I2_"
    else:
      print("WARNING: could not determine isospin")

    if "nonstrange" in self.direc:
      name += "S0_"
    elif "strange" in self.direc:
      name += "S1_"
    else:
      print("WARNING: could not determine strangeness")

    mom = self.direc[self.direc.find('mom'):].split('/')[0].split('_')
    if (len(mom) == 4):
      name += "P" + mom[1] + mom[2] + mom[3] + "_"
    else:
      print("WARNING: could not determine momentum")

    irrep = self.direc.rstrip('/').split('/')
    if (len(irrep)):
      irrep = irrep[len(irrep)-1]
      name += irrep
    else:
      print("WARNING: could not determine irrep")

    return name


  def __is_corr(self, filename):

    f = filename.split('/')
    f = f[len(f)-1]

    if f.find("corr_") == 0:
      return True
    else:
      return False
    
    '''
    result = subprocess.check_output([laph_query, filename]).strip()

    if result == "This is a LapH correlator file":
      return True
    else:
      return False
    '''
    
  def __is_vev(self, filename):

    f = filename.split('/')
    f = f[len(f)-1]
    
    if f.find("vev_") == 0:
      return True
    else:
      return False
    
    '''
    result = subprocess.check_output([laph_query, filename]).strip()

    if result == "This is a LapH VEV file":
      return True
    else:
      return False
    '''

  def __add_data_files(self, data_files, corrs):

    data_string = "VEVData"
    if corrs:
      data_string = "CorrelatorData"
      
    
    mc_tag = self.xml_root.find("Initialize").find("MCObservables")
    
    data_tag = mc_tag.find(data_string)
    if data_tag is None:
      data_tag = ET.SubElement(mc_tag, data_string)
      
    min_file = max_file = None
    file_stub = None
    for f in data_files:

      if file_stub is None:
        file_stub = f.split('.')[0]

      if min_file is None or max_file is None:
        min_file = max_file = int(f.split('.')[1])

      temp_file = int(f.split('.')[1])

      if temp_file < min_file:
        min_file = temp_file
      elif temp_file > max_file:
        max_file = temp_file

    file_list_info_tag = ET.SubElement(data_tag, "FileListInfo")
    file_name_stub_tag = ET.SubElement(file_list_info_tag, "FileNameStub")
    file_name_stub_tag.text = file_stub

    file_min_tag = ET.SubElement(file_list_info_tag, "MinFileNumber")
    file_min_tag.text = str(min_file)

    file_max_tag = ET.SubElement(file_list_info_tag, "MaxFileNumber")
    file_max_tag.text = str(max_file)


  def __set_operators(self, corr_files, vev_files, laph_query):

    for corr in corr_files:

      result = subprocess.check_output([laph_query, "-i", corr])
      result_xml = self.__laph_query_to_xml(result.decode())

      operator_1 = result_xml.find("CorrelatorInfo").find("Source").find("Operator")
      operator_2 = result_xml.find("CorrelatorInfo").find("Sink").find("Operator")

      self.ops.add(operators.Operator(operator_1))
      self.ops.add(operators.Operator(operator_2))

    for vev in vev_files:

      result = subprocess.check_output([laph_query, "-i", vev])
      result_xml = self.__laph_query_to_xml(result.decode())

      operator = result_xml.find("Operator")
      
      self.ops.add(operators.Operator(operator))

    
  def __laph_query_to_xml(self, laph_result):

    laph_result = laph_result[laph_result.find('<'):laph_result.rfind('>')+1]

    return ET.fromstring(laph_result)

