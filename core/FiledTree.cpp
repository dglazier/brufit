#include "FiledTree.h"
#include <TDirectory.h>

#include <utility>

namespace HS{
  namespace FIT{


    FiledTree::~FiledTree(){
      std::cout<<"FiledTree::~FiledTree()  tree name "<<fTree->GetName()<<" "<<fTree->GetEntries()<<" "<<fFile->GetName()<<endl;
      if(fMode==Mode_t::recreate||fMode==Mode_t::create||
	 fMode==Mode_t::update||fMode==Mode_t::copyfull||
	 fMode==Mode_t::copyempty){
	if(fFile){ 
	  fFile->cd();
	  if(fFile->IsWritable()){
	    Tree()->Write();
	  }
	  else{
	    std::cout<<"Warning deleting FiledTree but can't write to file!"<<endl;
	    std::cout<<"        this may mean you are letting ROOT cleanup "<<endl;
	    std::cout<<"        this object (i.e. after .q) rather than calling"<<endl;
	    std::cout<<"        its destructor somewhere in your code..."<<endl;
	    std::cout<<"       tree name "<<fTree->GetName()<<" "<<fFile->GetName()<<endl;
	  }
	}
      }
      fTree.reset();
    }
    filed_uptr FiledTree::Recreate(TString tname,TString fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(std::move(fname),"recreate");
      f->CreateTree(std::move(tname));
      f->SetMode(Mode_t::recreate);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::Read(const TString& tname,const TString& fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      auto file=TFile::Open(fname,"read");
      auto tree=dynamic_cast<TTree*>(file->Get(tname));
      f->SetFile(file);
      f->SetTree(tree);
      f->SetMode(Mode_t::read);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::Update(const TString& tname,const TString& fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      auto file=TFile::Open(fname,"update");
      auto tree=dynamic_cast<TTree*>(file->Get(tname));
      f->SetFile(file);
      f->SetTree(tree);
      f->SetMode(Mode_t::update);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::Create(TString tname,TString fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(std::move(fname),"create");
      f->CreateTree(std::move(tname));
      f->SetMode(Mode_t::create);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::CloneEmpty(TTree* tree,TString fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(std::move(fname),"recreate");
      f->SetTree(tree->CloneTree(0));
      f->SetMode(Mode_t::copyempty);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::CloneFull(TTree* tree,TString fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(std::move(fname),"recreate");
      f->SetTree(tree->CloneTree(-1,"fast"));
      f->SetMode(Mode_t::copyfull);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::CloneEmpty(const ttree_ptr& tree,TString fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(std::move(fname),"recreate");
      f->SetTree(tree->CloneTree(0));
      f->SetMode(Mode_t::copyempty);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::CloneFull(const ttree_ptr& tree,TString fname){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(std::move(fname),"recreate");
      f->SetTree(tree->CloneTree(-1,"fast"));
      f->SetMode(Mode_t::copyfull);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::RecreateCopyFull(const ttree_ptr& tree,const TString& fname,const TString& selection,const Long64_t nentries){
      //using copy tree allows use of tentrylists to filter
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(fname,"recreate");
      f->SetTree(tree->CopyTree(selection,"",nentries));
      f->SetMode(Mode_t::copyfull);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }
    filed_uptr FiledTree::RecreateCopyFull(TTree* tree,const TString& fname,const TString& selection,const Long64_t nentries ){
      auto saveDir=gDirectory;
      filed_uptr f{new FiledTree()};
      f->CreateFile(fname,"recreate");
      f->SetTree(tree->CopyTree(selection,"",nentries));
      f->SetMode(Mode_t::copyfull);
      f->SetTreeDirectory();
      saveDir->cd();
      return f;
    }

    
  } // namespace FIT
} // namespace HS
