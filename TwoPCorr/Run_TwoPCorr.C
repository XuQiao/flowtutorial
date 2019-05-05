void Run_TwoPCorr() {
    TwoPCorr* tp = new TwoPCorr();
    tp -> push_file("../Generator/vndata_1k_mult200.root");
    tp -> Init();
    tp -> ProcessEvents();
    tp -> End();
    tp -> Plot();
}
