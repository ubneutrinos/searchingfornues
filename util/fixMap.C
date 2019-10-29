/*
Usage:
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup root v6_18_04b -q e17:prof
root 'fixMap.C("file.root")'
*/
void fixMap(string filename){
    gInterpreter->GenerateDictionary("map<string,vector<double>>","map");
    
    string filenamemod = filename.substr(0,filename.find(".root"))+"_mapfixed.root";

    TFile file(filename.c_str());
    TFile newfile(filenamemod.c_str(),"recreate");

    TList* keylist = file.GetListOfKeys();
    TIter next(keylist);
    TKey* key;
    TObject* obj;

    while((key = (TKey*) next())){
        obj = key->ReadObj();
        if(obj->InheritsFrom("TDirectory")){
            TDirectory* newdir = newfile.mkdir(obj->GetName());
            TList* keylist2 = ((TDirectory*) obj)->GetListOfKeys();
            TIter next2(keylist2);
            TKey* key2;
            TObject* obj2;
                while((key2 = (TKey*) next2())){
                    obj2 = key2->ReadObj();
                    if(obj2->InheritsFrom("TTree")){
                        if((string)obj2->GetName() == "NeutrinoSelectionFilter"){
                            continue;
                        }
                        else{
                            TTree* newtree = ((TTree*)obj2)->CloneTree();
                            newdir->WriteObject(newtree,newtree->GetName());
                        }
                    }
                    else{
                        newdir->WriteObject(obj2,obj2->GetName());
                    }
                }
        }
        else{
            newfile.WriteObject(obj,obj->GetName());
        }
    }

    TTree* tree = nullptr;
    file.GetObject("nuselection/NeutrinoSelectionFilter",tree);

    tree->SetBranchStatus("weights",0);


    TDirectory* dir = newfile.GetDirectory("nuselection");
    dir->cd();

    TTree* newtree = tree->CloneTree(-1,"fast");

    tree->SetBranchStatus("weights",1);
    TBranch* b_weights = tree->GetBranch("weights");

    map<string,vector<double>>* weights = nullptr;    
    b_weights->SetAddress(&weights);
    
    TBranch* bnew = newtree->Branch("weights",weights);

    for(int i = 0; i<tree->GetEntries();++i){
        b_weights->GetEntry(i);
        bnew->Fill();
    }

    dir->WriteObject(newtree,newtree->GetName());
    gApplication->Terminate();
}
