# tRNA-charge-seq Architecture & Planning

**Purpose:** Reference architecture and design documents
**Last Updated:** 2026-02-10

---

## 📐 Package Architecture

- **[PACKAGE_ARCHITECTURE.md](PACKAGE_ARCHITECTURE.md)** - Package design reference
  - Dual interface (CLI + Python library)
  - minimap2-inspired design
  - Module organization
  - Nextflow integration

---

## 🎯 Current Development

**For active development plans, see:** `../.claude/CURRENT_PLAN.md`

**Key decisions:**
- Focus on charge quantification first
- Work with JSON.bz2 and CSV (native formats)
- Build lightweight alignment viewer (no IGV needed)
- 3-agent team for focused development

---

## 📝 Note

Original planning documents (20-week comprehensive plan, detailed implementation plans) have been superseded by the focused 8-week plan in `.claude/CURRENT_PLAN.md`.

The architecture design in this folder remains a useful reference for long-term package structure.
