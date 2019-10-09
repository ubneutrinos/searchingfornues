/* Usage:
root 'readMap.C("file.root")'
*/
void readMap(string filename){
   gInterpreter->GenerateDictionary("map<string,vector<double>>","map");

   TFile file(filename.c_str());

   TTree* tree = nullptr;
   file.GetObject("neutrinoselection_filt.root:/nuselection/NeutrinoSelectionFilter",tree);

   map<string,vector<double>>* weights = nullptr;
   TBranch* b_weights = tree->GetBranch("weights");
   
   b_weights->SetAddress(&weights);
   b_weights->GetEntry(0);

   for (auto& kv: *weights){
      for(auto& val: kv.second){
         cout<<kv.first<<" "<<val<<endl;
      }
   }
}
