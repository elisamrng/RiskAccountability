import xml
import xml.etree.ElementTree as ET


def clean(path_input_pnml, path_output_pnml):
    f = open(path_input_pnml, "r")
    x = open(path_output_pnml, "w")
    graphics = False
    line = f.readline()
    while line != '':
        if "<graphics>" not in line and graphics == False:
            if "skip" in line and "transition" in line:
                x.write(line) # transition
                x.write('<toolspecific tool="ProM" version="6.4" activity="$invisible$" localNodeID="45895f82-74eb-42c1-9ddf-e58c68f35d6a"/>')    
            else:
                x.write(line)
        else:
            if "<graphics>" in line:
                graphics = True
            if "</graphics>" in line:
                graphics = False
        line = f.readline()
    f.close()
    x.close()


# Renames the name of the transitions of a discovered model replacing
# white spaces with _
def rename_transitions(path_input_pnml,path_output_pnml):
    mytree = ET.parse(path_input_pnml)
    myroot = mytree.getroot()
    for el in myroot.iter():
        print(el)
    
    for t in myroot.iter('transition'):
        n = t.find("name").find("text")
        if (n.text).startswith("skip"):
            print("------------------")
            ET.SubElement(t,"prova")
            n.text = new_name
        t_id = t.get('id')
        t.set('id',n.text)

    mytree.write(path_output_pnml)




def make_align_txt_file(align, destination_path):
    f=open(destination_path,"w")
    list = align["alignment"]
    for item in list:
        if item[1]is None:
            f.write(item[0] + " " +"None\n")
        else:
            f.write(item[0] +" "+ item[1] + "\n")
    f.close()