#!/usr/bin/env python3
"""
Extract static HTML pages from Flask templates
Converts Jinja2 templates to standalone HTML files for nginx serving
"""

import sys
import re
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

def extract_static_pages(template_dir: Path, static_dir: Path, backend_dir: Path):
    """
    Extract static HTML pages from Flask templates
    
    Args:
        template_dir: Path to Flask templates directory
        static_dir: Output directory for static HTML files
        backend_dir: Path to backend directory (for static assets)
    """
    
    print("=" * 60)
    print("Extracting Static HTML Pages from Flask Templates")
    print("=" * 60)
    print(f"Template Directory: {template_dir}")
    print(f"Output Directory: {static_dir}")
    print("")
    
    # Create Jinja environment
    env = Environment(loader=FileSystemLoader(str(template_dir)))
    
    # Pages to convert to static HTML
    # Format: {template_filename: output_path}
    static_pages = {
        'index.html': 'index.html',
        'publications.html': 'publications/index.html',
        'research.html': 'research/index.html',
        'contact.html': 'contact/index.html',
    }
    
    # Process each template
    for template_name, output_path in static_pages.items():
        print(f"Processing: {template_name}")
        
        try:
            # Load template
            template = env.get_template(template_name)
            
            # Render template (with empty context)
            html_content = template.render()
            
            # Fix asset paths (Flask uses url_for, we need direct paths)
            html_content = fix_asset_paths(html_content)
            
            # Write to output file
            output_file = static_dir / output_path
            output_file.parent.mkdir(parents=True, exist_ok=True)
            output_file.write_text(html_content)
            
            print(f"  ✓ Created: {output_file}")
            print(f"  ✓ Size: {len(html_content)} bytes")
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
            continue
    
    print("")
    print("=" * 60)
    print("Static HTML extraction complete!")
    print("=" * 60)


def fix_asset_paths(html_content: str) -> str:
    """
    Fix asset paths in HTML for static serving
    
    Flask url_for() generates paths, but we need to ensure
    static assets point to /static/ correctly
    """
    
    # This function may not be needed if templates already use /static/ paths
    # But we keep it for potential future fixes
    
    # Example: Convert relative paths to absolute /static/ paths
    # html_content = html_content.replace('href="static/', 'href="/static/')
    # html_content = html_content.replace('src="static/', 'src="/static/')
    
    return html_content


def main():
    """Main entry point"""
    
    if len(sys.argv) < 3:
        print("Usage: python3 extract_static_pages.py <template_dir> <output_dir> [backend_dir]")
        print("")
        print("Example:")
        print("  python3 extract_static_pages.py \\")
        print("    /var/www/lab-website/backend/templates \\")
        print("    /var/www/lab-website/static-pages \\")
        print("    /var/www/lab-website/backend")
        sys.exit(1)
    
    template_dir = Path(sys.argv[1])
    static_dir = Path(sys.argv[2])
    backend_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else template_dir.parent
    
    # Validate paths
    if not template_dir.exists():
        print(f"Error: Template directory not found: {template_dir}")
        sys.exit(1)
    
    # Create output directory
    static_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract pages
    extract_static_pages(template_dir, static_dir, backend_dir)


if __name__ == '__main__':
    main()
