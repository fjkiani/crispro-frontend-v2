#!/usr/bin/env python3
"""
Attempt direct download of PTRC-HGSOC clinical data

Try various methods to download the CSV file directly.
"""

import requests
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

OUTPUT_DIR = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw" / "ptrc_hgsoc" / "clinical"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Known TCIA download patterns
TCIA_BASE = "https://www.cancerimagingarchive.net"
COLLECTION = "PTRC-HGSOC"

# Try various download URL patterns
DOWNLOAD_URLS = [
    # Direct CSV link patterns
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&SeriesInstanceUID=PTRC-HGSOC",
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&Collection={COLLECTION}",
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&CollectionName={COLLECTION}",
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&collection={COLLECTION.lower()}",
    # Clinical data specific
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&Collection={COLLECTION}&DataType=Clinical",
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&Collection={COLLECTION}&DataType=CSV",
    # Alternative patterns
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&Collection={COLLECTION}&format=csv",
    f"{TCIA_BASE}/nbia-search/services/SearchService?getfile&Collection={COLLECTION}&fileType=csv",
]

# Also try to get from collection page
COLLECTION_URL = f"{TCIA_BASE}/collection/ptrc-hgsoc/"


def try_download_url(url: str, filename: str) -> bool:
    """Try to download from a URL."""
    try:
        print(f"   Trying: {url[:80]}...")
        response = requests.get(url, timeout=30, allow_redirects=True)
        
        if response.status_code == 200:
            content_type = response.headers.get('Content-Type', '').lower()
            content = response.content
            
            # Check if it's CSV or text
            if 'csv' in content_type or 'text' in content_type or len(content) < 100000:
                # Try to decode as text
                try:
                    text = content.decode('utf-8')
                    # Check if it looks like CSV
                    if ',' in text or '\t' in text:
                        output_file = OUTPUT_DIR / filename
                        with open(output_file, 'wb') as f:
                            f.write(content)
                        print(f"   ✅ Downloaded: {output_file} ({len(content)} bytes)")
                        return True
                except:
                    pass
            
            # Save anyway if small file
            if len(content) < 1000000:  # Less than 1MB
                output_file = OUTPUT_DIR / filename
                with open(output_file, 'wb') as f:
                    f.write(content)
                print(f"   ✅ Downloaded: {output_file} ({len(content)} bytes)")
                return True
        
        return False
    
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return False


def extract_download_links_from_page():
    """Try to extract download links from the collection page."""
    print("\n📥 Extracting download links from collection page...")
    
    try:
        response = requests.get(COLLECTION_URL, timeout=30)
        if response.status_code == 200:
            content = response.text
            
            # Look for CSV links
            import re
            patterns = [
                r'href=["\']([^"\']*\.csv[^"\']*)["\']',
                r'href=["\']([^"\']*clinical[^"\']*[^"\']*)["\']',
                r'https://[^"\']*ptrc[^"\']*clinical[^"\']*',
            ]
            
            found_urls = set()
            for pattern in patterns:
                matches = re.findall(pattern, content, re.IGNORECASE)
                found_urls.update(matches)
            
            # Make URLs absolute
            absolute_urls = []
            for url in found_urls:
                if url.startswith('http'):
                    absolute_urls.append(url)
                elif url.startswith('/'):
                    absolute_urls.append(f"{TCIA_BASE}{url}")
                else:
                    absolute_urls.append(f"{COLLECTION_URL}{url}")
            
            return absolute_urls[:10]  # Limit to first 10
            
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return []


def main():
    """Try various methods to download the clinical data."""
    print("="*80)
    print("PTRC-HGSOC DIRECT DOWNLOAD ATTEMPT")
    print("="*80)
    print()
    
    # Try extracted URLs from page
    page_urls = extract_download_links_from_page()
    if page_urls:
        print(f"📥 Found {len(page_urls)} potential URLs from page")
        for i, url in enumerate(page_urls, 1):
            if try_download_url(url, f"clinical_data_from_page_{i}.csv"):
                print("   ✅ Success!")
                return
    
    # Try known URL patterns
    print("\n📥 Trying known URL patterns...")
    for i, url in enumerate(DOWNLOAD_URLS, 1):
        if try_download_url(url, f"clinical_data_pattern_{i}.csv"):
            print("   ✅ Success!")
            return
    
    print("\n" + "="*80)
    print("DOWNLOAD ATTEMPT FAILED")
    print("="*80)
    print()
    print("⚠️  Direct download not possible without authentication.")
    print("   TCIA requires login and manual download.")
    print()
    print("📋 Manual Download Required:")
    print(f"   1. Go to: {COLLECTION_URL}")
    print("   2. Login/register if needed")
    print("   3. Click 'Download' for Clinical data CSV")
    print(f"   4. Save to: {OUTPUT_DIR}")
    print()


if __name__ == "__main__":
    main()

