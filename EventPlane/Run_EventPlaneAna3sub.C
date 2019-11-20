void Run_EventPlaneAna3sub() {
	EventPlaneAna3sub* ep = new EventPlaneAna3sub();
	ep -> push_file("../Generator/vndata_1k_mult200.root");
	ep -> Init();
	ep -> set_nharm(2);
	ep -> ProcessEvents();
	ep -> End();
	ep -> Plot();
}
