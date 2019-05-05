void Run_QCumulant() {
	QCumulant* qcumu = new QCumulant();
	qcumu -> push_file("../Generator/vndata_1k_mult200.root");
	qcumu -> Init();
	qcumu -> set_nharm(2);
	qcumu -> ProcessEvents();
	qcumu -> End();
	qcumu -> Plot();
}
