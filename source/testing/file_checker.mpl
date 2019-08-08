with(StringTools):
interface(quiet=true):

fname1:="/latticeQCD/raid6/ahanlon/analysis/checks/32_B_2I1_S1_PSQ2_A1/32_B_2I1_S1_P0-11_A1_1_print_corr.log":
fname2:="/latticeQCD/raid6/ahanlon/analysis/checks/32_B_2I1_S1_PSQ2_A1/32_B_2I1_S1_P011_A1_1_print_corr.log":

for fname in [fname1,fname2] do

str:=readline(fname):
listen:=false:
donecorr:=false:
t35:=false:

while (str<>0) do

   if (Search("<Src>",str)>0) then
      src:=str:
      listen:=true:
      donecorr:=false:
      t35:=false:
   elif (Search("<Snk>",str)>0) then
      snk:=str:
      listen:=true:
      donecorr:=false:
      t35:=false:
   end if:
   if (Search("</PrintXML>",str)>0) then
      donecorr:=true:
      listen:=false:
      if (not t35) then
         printf("%s   %s\n",src,snk):
      end if:
      t35:=false:
   end if:
   if (listen) then
      if (Search("<TimeSeparation>35",str)>0) then
         t35:=true:
      end if:
   end if:
   str:=readline(fname):
   end do:

printf("******************\n"):
end do:
