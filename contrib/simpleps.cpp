/* Simple postscript library - Matteo Cicuttin (c) 2016 */

#include "simpleps.h"

namespace sps {

path::path()
{}

void
path::emit_code(std::ostream& os)
{
    if (m_points.size() < 2)
        return;

    auto firstpoint = *m_points.begin();

    os << "newpath" << std::endl;
    os << firstpoint.x() << " " << firstpoint.y() << " moveto" << std::endl;
    for (auto itor = std::next(m_points.begin()); itor != m_points.end(); itor++)
    {
        auto curpoint = *itor;
        os << curpoint.x() << " " << curpoint.y() << " lineto" << std::endl;
    }
    os << "stroke" << std::endl;
}

text::text()
{}

void
text::emit_code(std::ostream& os)
{
    os << "/Times-Roman findfont" << std::endl;
    os << "5 scalefont" << std::endl;
    os << "setfont" << std::endl;
    os << "newpath" << std::endl;
    os << m_position.x() << " " << m_position.y() << " moveto" << std::endl;
    os << "(" << m_text << ") show" << std::endl;
}

simpleps::simpleps()
{}

void
simpleps::add(postscript_object *po)
{
    m_ps_objs.push_back(po);
}

bool
simpleps::write(const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        return false;

    for (auto& po : m_ps_objs)
        po->emit_code(ofs);

    ofs.close();
    return true;
}

void
simpleps::output(std::ostream& os)
{
    for (auto& po : m_ps_objs)
        po->emit_code(os);
}

simpleps::~simpleps()
{
    for (auto& po : m_ps_objs)
        delete po;
}

} // namespace simpleps
