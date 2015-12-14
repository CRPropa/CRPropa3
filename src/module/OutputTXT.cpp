#include "crpropa/module/OutputTXT.h"
#include "crpropa/Units.h"

#include <stdio.h>
#include <stdexcept>

#ifdef CRPROPA_HAVE_ZLIB
#include <ozstream.hpp>
#include <kiss/string.h>
#endif

namespace crpropa {

TextOutput::TextOutput() :
		lengthScale(Mpc), energyScale(EeV), out(&std::cout), oneDimensional(
				false) {
	enableAll();
}

TextOutput::TextOutput(std::ostream &out) :
		lengthScale(Mpc), energyScale(EeV), out(&out), oneDimensional(false) {
	enableAll();
}

TextOutput::TextOutput(const std::string &filename) :
		lengthScale(Mpc), energyScale(EeV), outfile(filename.c_str(),
				std::ios::binary), out(&outfile), oneDimensional(false), filename(
				filename) {
	enableAll();
	if (kiss::ends_with(filename, ".gz"))
		gzip();
}

TextOutput::TextOutput(const std::string &filename,
		const std::string &outputtype) :
		lengthScale(Mpc), energyScale(EeV), outfile(filename.c_str(),
				std::ios::binary), out(&outfile), oneDimensional(false), filename(
				filename) {
	if (outputtype == "1D trajectories") {
		// X, ID, E
		set(CurrentPositionColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set1D(true);
	} else if (outputtype == "1D events") {
		// ID, E, D, ID0, E0
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(TrajectoryLengthColumn, true);
		set(SourceIdColumn, true);
		set(SourceEnergyColumn, true);
		set1D(true);
	} else if (outputtype == "3D trajectories") {
		// D, ID, E, X, Y, Z, Px, Py, Pz
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(CurrentPositionColumn, true);
		set(CurrentDirectionColumn, true);
	} else if (outputtype == "3D events") {
		// D, ID, ID0, E, E0, X, Y, Z, X0, Y0, Z0, Px, Py, Pz, P0x, P0y, P0z,z
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(SourceIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(SourceEnergyColumn, true);
		set(CurrentPositionColumn, true);
		set(SourcePositionColumn, true);
		set(CurrentDirectionColumn, true);
		set(SourcePositionColumn, true);
	} else {
		throw std::runtime_error(
				"TextOutput: Outputtype must be one of '1D trajectories', '1D events', '3D trajectories', '3D events'");
	}
	if (kiss::ends_with(filename, ".gz"))
		gzip();
}

void TextOutput::setEnergyScale(double scale) {
	energyScale = scale;
	printHeader();
}

void TextOutput::setLengthScale(double scale) {
	lengthScale = scale;
	printHeader();
}

void TextOutput::set1D(bool value) {
	oneDimensional = value;
	printHeader();
}

void TextOutput::enable(OutputColumn field) {
	fields.set(field, true);
	printHeader();
}

void TextOutput::disable(OutputColumn field) {
	fields.set(field, false);
	printHeader();
}

void TextOutput::set(OutputColumn field, bool value) {
	fields.set(field, value);
	printHeader();
}

void TextOutput::enableAll() {
	fields.set();
	printHeader();
}

void TextOutput::disableAll() {
	fields.reset();
	printHeader();
}

void TextOutput::resetFile() {
	// clear output file
	outfile.close();
	outfile.open(filename.c_str(), std::ios::binary);
}

void TextOutput::printHeader() {
	resetFile();

	*out << "#";
	if (fields.test(TrajectoryLengthColumn))
		*out << "\tD";
	if (fields.test(RedshiftColumn))
		*out << "\tz";
	if (fields.test(CurrentIdColumn))
		*out << "\tID";
	if (fields.test(CurrentEnergyColumn))
		*out << "\tE";
	if (fields.test(CurrentPositionColumn)) {
		if (oneDimensional)
			*out << "\tX";
		else
			*out << "\tX\tY\tZ";
	}
	if (fields.test(CurrentDirectionColumn))
		if (not oneDimensional)
			*out << "\tPx\tPy\tPz";

	if (fields.test(SourceIdColumn))
		*out << "\tID0";
	if (fields.test(SourceEnergyColumn))
		*out << "\tE0";
	if (fields.test(SourcePositionColumn)) {
		if (oneDimensional)
			*out << "\tX0";
		else
			*out << "\tX0\tY0\tZ0";
	}
	if (fields.test(SourceDirectionColumn))
		if (not oneDimensional)
			*out << "\tP0x\tP0y\tP0z";

	if (fields.test(CreatedIdColumn))
		*out << "\tID1";
	if (fields.test(CreatedEnergyColumn))
		*out << "\tE1";
	if (fields.test(CreatedPositionColumn)) {
		if (oneDimensional)
			*out << "\tX1";
		else
			*out << "\tX1\tY1\tZ1";
	}
	if (fields.test(CreatedDirectionColumn))
		if (not oneDimensional)
			*out << "\tP1x\tP1y\tP1z";

	*out << "\n#\n";
	if (fields.test(TrajectoryLengthColumn))
		*out << "# D             Trajectory length [" << lengthScale / Mpc
				<< " Mpc]\n";
	if (fields.test(RedshiftColumn))
		*out << "# z             Redshift\n";
	if (fields.test(CurrentIdColumn) || fields.test(CreatedIdColumn)
			|| fields.test(SourceIdColumn))
		*out << "# ID/ID0/ID1    Particle type (PDG MC numbering scheme)\n";
	if (fields.test(CurrentEnergyColumn) || fields.test(CreatedEnergyColumn)
			|| fields.test(SourceEnergyColumn))
		*out << "# E/E0/E1       Energy [" << energyScale / EeV << " EeV]\n";
	if (fields.test(CurrentPositionColumn) || fields.test(CreatedPositionColumn)
			|| fields.test(SourcePositionColumn))
		*out << "# X/X0/X1...    Position [" << lengthScale / Mpc << " Mpc]\n";
	if (fields.test(CurrentDirectionColumn)
			|| fields.test(CreatedDirectionColumn)
			|| fields.test(SourceDirectionColumn))
		*out << "# Px/P0x/P1x... Heading (unit vector of momentum)\n";
	*out
			<< "# no index = current, 0 = at source, 1 = at point of creation\n#\n";
}

void TextOutput::process(Candidate *c) const {
	if (fields.none())
		return;

	char buffer[1024];
	size_t p = 0;

	if (fields.test(TrajectoryLengthColumn))
		p += sprintf(buffer + p, "%8.5f\t",
				c->getTrajectoryLength() / lengthScale);
	if (fields.test(RedshiftColumn))
		p += sprintf(buffer + p, "%1.5f\t", c->getRedshift());
	if (fields.test(CurrentIdColumn))
		p += sprintf(buffer + p, "%10i\t", c->current.getId());
	if (fields.test(CurrentEnergyColumn))
		p += sprintf(buffer + p, "%8.5f\t",
				c->current.getEnergy() / energyScale);
	if (fields.test(CurrentPositionColumn)) {
		if (oneDimensional) {
			p += sprintf(buffer + p, "%8.5f\t",
					c->current.getPosition().x / lengthScale);
		} else {
			const Vector3d pos = c->current.getPosition() / lengthScale;
			p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", pos.x, pos.y,
					pos.z);
		}
	}
	if (fields.test(CurrentDirectionColumn)) {
		if (not oneDimensional) {
			const Vector3d pos = c->current.getDirection();
			p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", pos.x, pos.y,
					pos.z);
		}
	}

	if (fields.test(SourceIdColumn))
		p += sprintf(buffer + p, "%10i\t", c->source.getId());
	if (fields.test(SourceEnergyColumn))
		p += sprintf(buffer + p, "%8.5f\t",
				c->source.getEnergy() / energyScale);
	if (fields.test(SourcePositionColumn)) {
		if (oneDimensional) {
			p += sprintf(buffer + p, "%8.5f\t",
					c->source.getPosition().x / lengthScale);
		} else {
			const Vector3d pos = c->source.getPosition() / lengthScale;
			p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", pos.x, pos.y,
					pos.z);
		}
	}
	if (fields.test(SourceDirectionColumn)) {
		if (not oneDimensional) {
			const Vector3d pos = c->source.getDirection();
			p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", pos.x, pos.y,
					pos.z);
		}

	}

	if (fields.test(CreatedIdColumn))
		p += sprintf(buffer + p, "%10i\t", c->created.getId());
	if (fields.test(CreatedEnergyColumn))
		p += sprintf(buffer + p, "%8.5f\t",
				c->created.getEnergy() / energyScale);
	if (fields.test(CreatedPositionColumn)) {
		if (oneDimensional) {
			p += sprintf(buffer + p, "%8.5f\t",
					c->created.getPosition().x / lengthScale);
		} else {
			const Vector3d pos = c->created.getPosition() / lengthScale;
			p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", pos.x, pos.y,
					pos.z);
		}
	}
	if (fields.test(CreatedDirectionColumn)) {
		if (not oneDimensional) {
			const Vector3d pos = c->created.getDirection();
			p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", pos.x, pos.y,
					pos.z);
		}
	}

	buffer[p - 1] = '\n';

#pragma omp critical
	{
		out->write(buffer, p);
	}

}

std::string TextOutput::getDescription() const {
	return "TextOutput";
}


void TextOutput::endRun() {
#ifdef CRPROPA_HAVE_ZLIB
	zstream::ogzstream *zs = dynamic_cast<zstream::ogzstream *>(out);
	if (zs) {
		zs->close();
	}
#endif    
}

TextOutput::~TextOutput() {
#ifdef CRPROPA_HAVE_ZLIB
	zstream::ogzstream *zs = dynamic_cast<zstream::ogzstream *>(out);
	if (zs) {
		delete zs;
	}
#endif
	outfile.close();
}

void TextOutput::gzip() {
#ifdef CRPROPA_HAVE_ZLIB
	out = new zstream::ogzstream(*out);
#else
	throw std::runtime_error("CRPropa was build without Zlib compression!");
#endif
}

TrajectoryOutput::TrajectoryOutput(std::string name) {
	std::cout
			<< "TrajectoryOutput is deprecated. Use TextOutput(..., '3D trajectories') instead."
			<< std::endl;
	setDescription("Trajectory output");
	fout.open(name.c_str());
	fout << "# D\tID\tE\tX\tY\tZ\tPx\tPy\tPz\n" << "#\n"
			<< "# D           Trajectory length\n"
			<< "# ID          Particle type (PDG MC numbering scheme)\n"
			<< "# E           Energy [EeV]\n"
			<< "# X, Y, Z     Position [Mpc]\n"
			<< "# Px, Py, Pz  Heading (unit vector of momentum)\n" << "#\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	fout.close();
}

void TrajectoryOutput::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%8.3f\t", c->getTrajectoryLength() / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	Vector3d pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%8.4f\t%8.4f\t%8.4f\t", pos.x, pos.y, pos.z);
	const Vector3d &dir = c->current.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\n", dir.x, dir.y, dir.z);

#pragma omp critical
	fout.write(buffer, p);
}

void TrajectoryOutput::endRun() {
	fout.flush();
}

ConditionalOutput::ConditionalOutput(std::string fname, std::string cond) :
		condition(cond) {
	std::cout
			<< "ConditionalOutput is deprecated. Use TextOutput(..., '3D events') instead."
			<< std::endl;
	setDescription(
			"Conditional output, condition: " + cond + ", filename: " + fname);
	fout.open(fname.c_str());
	fout
			<< "# D\tID\tID0\tE\tE0\tX\tY\tZ\tX0\tY0\tZ0\tPx\tPy\tPz\tP0x\tP0y\tP0z\tz\n"
			<< "#\n" << "# D           Trajectory length [Mpc]\n"
			<< "# ID          Particle type (PDG MC numbering scheme)\n"
			<< "# E           Energy [EeV]\n"
			<< "# X, Y, Z     Position [Mpc]\n"
			<< "# Px, Py, Pz  Heading (unit vector of momentum)\n"
			<< "# z           Current redshift\n"
			<< "# Initial state: ID0, E0, ...\n" << "#\n";
}

ConditionalOutput::~ConditionalOutput() {
	fout.close();
}

void ConditionalOutput::process(Candidate *c) const {
	if (not (c->hasProperty(condition)))
		return;

	c->removeProperty(condition);

	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%8.3f\t", c->getTrajectoryLength() / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%10i\t", c->source.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	p += sprintf(buffer + p, "%8.4f\t", c->source.getEnergy() / EeV);
	Vector3d pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", pos.x, pos.y, pos.z);
	Vector3d ipos = c->source.getPosition() / Mpc;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", ipos.x, ipos.y, ipos.z);
	Vector3d dir = c->current.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", dir.x, dir.y, dir.z);
	Vector3d idir = c->source.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", idir.x, idir.y, idir.z);
	p += sprintf(buffer + p, "%1.3f\n", c->getRedshift());

#pragma omp critical
	fout.write(buffer, p);
}

void ConditionalOutput::endRun() {
	fout.flush();
}

TrajectoryOutput1D::TrajectoryOutput1D(std::string filename) {
	std::cout
			<< "TrajectoryOutput1D is deprecated. Use TextOutput(..., '1D trajectories') instead."
			<< std::endl;
	setDescription("TrajectoryOutput, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#X\tID\tE\n" << "#\n" << "# X  Position [Mpc]\n"
			<< "# ID Particle type\n" << "# E  Energy [EeV]\n";
}

TrajectoryOutput1D::~TrajectoryOutput1D() {
	fout.close();
}

void TrajectoryOutput1D::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;
	p += sprintf(buffer + p, "%8.4f\t", c->current.getPosition().x / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\n", c->current.getEnergy() / EeV);

#pragma omp critical
	fout.write(buffer, p);
}

void TrajectoryOutput1D::endRun() {
	fout.flush();
}

EventOutput1D::EventOutput1D(std::string filename) {
	std::cout
			<< "EventOutput1D is deprecated. Use TextOutput(..., '1D events') instead."
			<< std::endl;
	setDescription("Conditional output, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#ID\tE\tD\tID0\tE0\n" << "#\n" << "# ID  Particle type\n"
			<< "# E   Energy [EeV]\n"
			<< "# D   Comoving source distance [Mpc]\n"
			<< "# ID0 Initial particle type\n"
			<< "# E0  Initial energy [EeV]\n";
}

EventOutput1D::~EventOutput1D() {
	fout.close();
}

void EventOutput1D::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	p += sprintf(buffer + p, "%9.4f\t", c->source.getPosition().x / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->source.getId());
	p += sprintf(buffer + p, "%8.4f\n", c->source.getEnergy() / EeV);

#pragma omp critical
	fout.write(buffer, p);
}

void EventOutput1D::endRun() {
	fout.flush();
}

} // namespace crpropa
