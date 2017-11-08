#include "crpropa/module/TextOutput.h"
#include "crpropa/module/ParticleCollector.h"
#include "crpropa/Units.h"

#include <stdio.h>
#include <stdexcept>
#include <iostream>
#include <kiss/string.h>

#ifdef CRPROPA_HAVE_ZLIB
#include <izstream.hpp>
#include <ozstream.hpp>
#endif

namespace crpropa {

TextOutput::TextOutput() : Output(), out(&std::cout) {
}

TextOutput::TextOutput(OutputType outputtype) : Output(outputtype), out(&std::cout) {
}

TextOutput::TextOutput(std::ostream &out) : Output(), out(&out) {

}

TextOutput::TextOutput(std::ostream &out,
		OutputType outputtype) : Output(outputtype), out(&out) {
}

TextOutput::TextOutput(const std::string &filename) :  Output(), outfile(filename.c_str(),
				std::ios::binary), out(&outfile),  filename(
				filename) {
	if (kiss::ends_with(filename, ".gz"))
		gzip();
}

TextOutput::TextOutput(const std::string &filename,
				OutputType outputtype) : Output(outputtype), outfile(filename.c_str(),
				std::ios::binary), out(&outfile), filename(
				filename) {
	if (kiss::ends_with(filename, ".gz"))
		gzip();
}


void TextOutput::printHeader() const {
	*out << "#";
	if (fields.test(TrajectoryLengthColumn))
		*out << "\tD";
	if (fields.test(RedshiftColumn))
		*out << "\tz";
	if (fields.test(SerialNumberColumn))
		*out << "\tSN";
	if (fields.test(CurrentIdColumn))
		*out << "\tID";
	if (fields.test(CurrentEnergyColumn))
		*out << "\tE";
	if (fields.test(CurrentPositionColumn) && oneDimensional)
		*out << "\tX";
	if (fields.test(CurrentPositionColumn) && not oneDimensional)
		*out << "\tX\tY\tZ";
	if (fields.test(CurrentDirectionColumn) && not oneDimensional)
		*out << "\tPx\tPy\tPz";
	if (fields.test(SerialNumberColumn))
		*out << "\tSN0";
	if (fields.test(SourceIdColumn))
		*out << "\tID0";
	if (fields.test(SourceEnergyColumn))
		*out << "\tE0";
	if (fields.test(SourcePositionColumn) && oneDimensional) 
		*out << "\tX0";
	if (fields.test(SourcePositionColumn) && not oneDimensional)
		*out << "\tX0\tY0\tZ0";
	if (fields.test(SourceDirectionColumn) && not oneDimensional)
		*out << "\tP0x\tP0y\tP0z";
	if (fields.test(SerialNumberColumn))
		*out << "\tSN1";
	if (fields.test(CreatedIdColumn))
		*out << "\tID1";
	if (fields.test(CreatedEnergyColumn))
		*out << "\tE1";
	if (fields.test(CreatedPositionColumn) && oneDimensional)
		*out << "\tX1";
	if (fields.test(CreatedPositionColumn) && not oneDimensional)
		*out << "\tX1\tY1\tZ1";
	if (fields.test(CreatedDirectionColumn) && not oneDimensional)
		*out << "\tP1x\tP1y\tP1z";
	if (fields.test(WeightColumn))
		*out << "\tW";
	for(std::vector<Property>::const_iterator iter = properties.begin();
			iter != properties.end(); ++iter)
	{
		*out << "\t" << (*iter).name;
	}

	*out << "\n#\n";
	if (fields.test(TrajectoryLengthColumn))
		*out << "# D             Trajectory length [" << lengthScale / Mpc
				<< " Mpc]\n";
	if (fields.test(RedshiftColumn))
		*out << "# z             Redshift\n";
	if (fields.test(SerialNumberColumn))
		*out << "# SN/SN0/SN1    Serial number. Unique (within this run) id of the particle.\n";
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
	if (fields.test(WeightColumn))
		*out << "# W             Weights" << " \n";
	for(std::vector<Property>::const_iterator iter = properties.begin();
			iter != properties.end(); ++iter)
	{
			*out << "# " << (*iter).name << " " << (*iter).comment << "\n";
	}

	*out << "# no index = current, 0 = at source, 1 = at point of creation\n#\n";
}

void TextOutput::process(Candidate *c) const {
	if (fields.none() && properties.empty())
		return;

	char buffer[1024];
	size_t p = 0;

	std::locale old_locale = std::locale::global(std::locale::classic());

	if (fields.test(TrajectoryLengthColumn))
		p += sprintf(buffer + p, "%8.5E\t",
				c->getTrajectoryLength() / lengthScale);
	if (fields.test(RedshiftColumn))
		p += sprintf(buffer + p, "%1.5E\t", c->getRedshift());

	if (fields.test(SerialNumberColumn))
		p += sprintf(buffer + p, "%10lu\t",
				c->getSerialNumber());
	if (fields.test(CurrentIdColumn))
		p += sprintf(buffer + p, "%10i\t", c->current.getId());
	if (fields.test(CurrentEnergyColumn))
		p += sprintf(buffer + p, "%8.5E\t",
				c->current.getEnergy() / energyScale);
	if (fields.test(CurrentPositionColumn)) {
		if (oneDimensional) {
			p += sprintf(buffer + p, "%8.5E\t",
					c->current.getPosition().x / lengthScale);
		} else {
			const Vector3d pos = c->current.getPosition() / lengthScale;
			p += sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos.x, pos.y,
					pos.z);
		}
	}
	if (fields.test(CurrentDirectionColumn)) {
		if (not oneDimensional) {
			const Vector3d pos = c->current.getDirection();
			p += sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos.x, pos.y,
					pos.z);
		}
	}

	if (fields.test(SerialNumberColumn))
		p += sprintf(buffer + p, "%10lu\t", c->getSourceSerialNumber());
	if (fields.test(SourceIdColumn))
		p += sprintf(buffer + p, "%10i\t", c->source.getId());
	if (fields.test(SourceEnergyColumn))
		p += sprintf(buffer + p, "%8.5E\t",
				c->source.getEnergy() / energyScale);
	if (fields.test(SourcePositionColumn)) {
		if (oneDimensional) {
			p += sprintf(buffer + p, "%8.5E\t",
					c->source.getPosition().x / lengthScale);
		} else {
			const Vector3d pos = c->source.getPosition() / lengthScale;
			p += sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos.x, pos.y,
					pos.z);
		}
	}
	if (fields.test(SourceDirectionColumn)) {
		if (not oneDimensional) {
			const Vector3d pos = c->source.getDirection();
			p += sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos.x, pos.y,
					pos.z);
		}

	}

	if (fields.test(SerialNumberColumn))
		p += sprintf(buffer + p, "%10lu\t",
				c->getCreatedSerialNumber());
	if (fields.test(CreatedIdColumn))
		p += sprintf(buffer + p, "%10i\t", c->created.getId());
	if (fields.test(CreatedEnergyColumn))
		p += sprintf(buffer + p, "%8.5E\t",
				c->created.getEnergy() / energyScale);
	if (fields.test(CreatedPositionColumn)) {
		if (oneDimensional) {
			p += sprintf(buffer + p, "%8.5E\t",
					c->created.getPosition().x / lengthScale);
		} else {
			const Vector3d pos = c->created.getPosition() / lengthScale;
			p += sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos.x, pos.y,
					pos.z);
		}
	}
	if (fields.test(CreatedDirectionColumn)) {
		if (not oneDimensional) {
			const Vector3d pos = c->created.getDirection();
			p += sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos.x, pos.y,
					pos.z);
		}
	}
	if (fields.test(WeightColumn)) {
		p += sprintf(buffer + p, "%8.5E\t", c->getWeight());
	}

	for(std::vector<Output::Property>::const_iterator iter = properties.begin();
			iter != properties.end(); ++iter)
	{
		  Variant v;
			if (c->hasProperty((*iter).name))
			{
				v = c->getProperty((*iter).name);
			}
			else
			{
				v = (*iter).defaultValue;
			}
			p += sprintf(buffer + p, v.toString().c_str());
			p += sprintf(buffer + p, "\t");
	}
	buffer[p - 1] = '\n';

	std::locale::global(old_locale);

#pragma omp critical
	{
		if (count == 0)
			printHeader();
		Output::process(c);
		out->write(buffer, p);
	}

}

void TextOutput::load(const std::string &filename, ParticleCollector *collector){

	std::string line;
	std::istream *in;
	std::ifstream infile(filename.c_str());
	
	double lengthScale = Mpc; // default Mpc
	double energyScale = EeV; // default EeV

	if (!infile.good())
		throw std::runtime_error("crpropa::TextOutput: could not open file " + filename);
	in = &infile;
	
	if (kiss::ends_with(filename, ".gz")){
#ifdef CRPROPA_HAVE_ZLIB
		in = new zstream::igzstream(*in);
#else
		throw std::runtime_error("CRPropa was build without Zlib compression!");
#endif
	}

	while (std::getline(*in,line)) {
		std::stringstream stream(line);
		if (stream.peek() == '#')
			continue;

		ref_ptr<Candidate> c = new Candidate(); 
		double val_d; int val_i;
		double x, y, z;
		stream >> val_d;
		c->setTrajectoryLength(val_d*lengthScale); // D
		stream >> val_d;
		c->setRedshift(val_d); // z
		stream >> val_i;
		c->setSerialNumber(val_i); // SN
		stream >> val_i;
		c->current.setId(val_i); // ID
		stream >> val_d;
		c->current.setEnergy(val_d*energyScale); // E
		stream >> x >> y >> z;
		c->current.setPosition(Vector3d(x, y, z)*lengthScale); // X, Y, Z
		stream >> x >> y >> z;
		c->current.setDirection(Vector3d(x, y, z)*lengthScale); // Px, Py, Pz
		stream >> val_i; // SN0 (TODO: Reconstruct the parent-child relationship)
		stream >> val_i;
		c->source.setId(val_i); // ID0
		stream >> val_d;
		c->source.setEnergy(val_d*energyScale);	// E0
		stream >> x >> y >> z;
		c->source.setPosition(Vector3d(x, y, z)*lengthScale); // X0, Y0, Z0
		stream >> x >> y >> z;
		c->source.setDirection(Vector3d(x, y, z)*lengthScale); // P0x, P0y, P0z
		stream >> val_i; // SN1
		stream >> val_i;
		c->created.setId(val_i); // ID1
		stream >> val_d;
		c->created.setEnergy(val_d*energyScale); // E1
		stream >> x >> y >> z;
		c->created.setPosition(Vector3d(x, y, z)*lengthScale); // X1, Y1, Z1
		stream >> x >> y >> z;
		c->created.setDirection(Vector3d(x, y, z)*lengthScale); // P1x, P1y, P1z
		stream >> val_d;
		c->setWeight(val_d); // W

		collector->process(c);
	}
	infile.close();
}

std::string TextOutput::getDescription() const {
	return "TextOutput";
}

void TextOutput::close() {
#ifdef CRPROPA_HAVE_ZLIB
	zstream::ogzstream *zs = dynamic_cast<zstream::ogzstream *>(out);
	if (zs) {
		zs->close();
		delete out;
		out = 0;
	}
#endif
	outfile.flush();
}

TextOutput::~TextOutput() {
	close();
}

void TextOutput::gzip() {
#ifdef CRPROPA_HAVE_ZLIB
	out = new zstream::ogzstream(*out);
#else
	throw std::runtime_error("CRPropa was build without Zlib compression!");
#endif
}

} // namespace crpropa
