// This file contains the documentation for the DANA framework
// It is written in this way so Doxygen can find it.

/**
 * @defgroup DANA the DANA framework
 * @{
*/
		class DContainer{};
		class DEvent{};
		class DEventLoop[];
		class DEventProcessor{};
		class DEventSource{};
		class DEventSourceET{};
		class DEventSourceFile{};
		class DFactory{};
/**

	\page DANA The DANA framework

	\section Introduction
	
	The DANA (Hall-<b>D</b> <b>Ana</b>lysis) framework is a framework
	which the Hall-D reconstruction/analysis software is built around.
	The framework implements the HDDM (<b>H</b>all-<b>D</b> <b>D</b>ata
	<b>M</b>odel) I/O package for reading in HDDM formatted files. The
	framework is designed however so that other input sources can be used.
	<p></p>
	The framework itself is designed to work with event based data.
	A program typically creates a single DEventLoop object which
	implements the main event loop for the program. User routines
	are invoked as <i>callbacks</i> during the event loop. Reconstruction
	routines are implemented as <i>factories</i> which can be accessed
	from the callback routines. By implementing the reconstruction as
	factories, reconstruction code is executed only when it is needed.
	Having partial reconstruction for some events allows for efficient
	filtering of events as well as more efficient "special purpose"
	programs which look only at specific detector systems.
	

	@}
*/
