void Run_EventPlaneAna() {
	EventPlaneAna* ep = new EventPlaneAna();
	ep -> push_file("../Generator/vndata_1k_mult200.root");
	ep -> Init();
	ep -> set_nharm(2);
	ep -> ProcessEvents();
	ep -> End();
	ep -> Plot();
}
