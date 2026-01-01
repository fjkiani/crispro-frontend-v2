#!/usr/bin/env python3
"""Project Data Sphere API Client for SAS CAS API access."""

import os
from typing import List, Dict, Optional

try:
    import swat
    SWAT_AVAILABLE = True
except ImportError:
    SWAT_AVAILABLE = False


class ProjectDataSphereClient:
    def __init__(self, cas_url=None, ssl_cert_path=None, port=443):
        if not SWAT_AVAILABLE:
            raise ImportError("swat library required. Install with: pip install swat")
        
        self.cas_url = cas_url
        self.ssl_cert_path = ssl_cert_path
        self.port = port
        self.conn = None
        
        if ssl_cert_path and os.path.exists(ssl_cert_path):
            os.environ["CAS_CLIENT_SSL_CA_LIST"] = ssl_cert_path
            print(f"✅ SSL certificate set: {ssl_cert_path}")
    
    def connect(self, username=None, password=None):
        if not self.cas_url:
            print("❌ CAS URL not set")
            return False
        
        try:
            print(f"🔌 Connecting to: {self.cas_url}:{self.port}")
            if username and password:
                self.conn = swat.CAS(self.cas_url, self.port, username=username, password=password)
            else:
                self.conn = swat.CAS(self.cas_url, self.port)
            
            status = self.conn.serverstatus()
            print(f"✅ Connected! Status: {status}")
            return True
        except Exception as e:
            print(f"❌ Connection failed: {e}")
            return False
    
    def list_caslibs(self):
        if not self.conn:
            return []
        try:
            result = self.conn.table.caslibInfo()
            caslibs = []
            if hasattr(result, 'CASLibInfo'):
                for caslib in result.CASLibInfo:
                    caslibs.append({
                        'name': caslib.get('Name', ''),
                        'description': caslib.get('Description', ''),
                        'path': caslib.get('Path', '')
                    })
            return caslibs
        except Exception as e:
            print(f"❌ Error: {e}")
            return []
    
    def list_files_in_caslib(self, caslib_name, path=None):
        if not self.conn:
            return []
        try:
            params = {'allFiles': True, 'caslib': caslib_name}
            if path:
                params['path'] = path
            result = self.conn.table.fileInfo(**params)
            files = []
            if hasattr(result, 'FileInfo'):
                for f in result.FileInfo:
                    files.append({
                        'name': f.get('Name', ''),
                        'path': f.get('Path', ''),
                        'size': f.get('Size', 0)
                    })
            return files
        except Exception as e:
            print(f"❌ Error: {e}")
            return []
    
    
    
    def search_for_ovarian_cancer_data(self):
        if not self.conn:
            return []
        ovarian_datasets = []
        caslibs = self.list_caslibs()
        print(f"\n🔍 Searching {len(caslibs)} caslibs...")
        keywords = ['ovarian', 'ovary', 'ov', 'ca125', 'ca-125', 'platinum', 'pfi', 'serous', 'hgsc']
        for caslib in caslibs:
            caslib_name = caslib.get('name', '')
            files = self.list_files_in_caslib(caslib_name)
            for f in files:
                name_lower = f.get('name', '').lower()
                for keyword in keywords:
                    if keyword in name_lower:
                        ovarian_datasets.append({
                            'caslib': caslib_name,
                            'file_name': f.get('name', ''),
                            'matched_keyword': keyword
                        })
                        print(f"   ✅ Found: {f.get('name', '')} ({keyword})")
                        break
        return ovarian_datasets
    
    def disconnect(self):
        if self.conn:
            try:
                self.conn.close()
                print("✅ Disconnected")
            except:
                pass
            self.conn = None


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cas-url", default="https://mpmprodvdmml.ondemand.sas.com/cas-shared-default-http/")
    parser.add_argument("--ssl-cert", default="/Users/fahadkiani/Desktop/development/crispr-assistant-main/data/certs/trustedcerts.pem")
    parser.add_argument("--username", default="mpm0fxk2", help="SAS username")
    parser.add_argument("--password", help="PDS password (same as SAS password)")
    args = parser.parse_args()
    
    print("=" * 80)
    print("Project Data Sphere Connection Test")
    print("=" * 80)
    
    client = ProjectDataSphereClient(cas_url=args.cas_url, ssl_cert_path=args.ssl_cert)
    
    if client.connect(username=args.username, password=args.password):
        print("\nAvailable Caslibs:")
        caslibs = client.list_caslibs()
        for c in caslibs[:20]:
            print(f"  - {c.get('name', 'N/A')}")
        
        print("\nSearching for Ovarian Cancer Data:")
        ovarian_data = client.search_for_ovarian_cancer_data()
        
        if ovarian_data:
            print(f"\n✅ Found {len(ovarian_data)} datasets:")
            for d in ovarian_data:
                print(f"   - {d['caslib']}/{d['file_name']}")
        else:
            print("\n⚠️  No ovarian cancer datasets found")
        
        client.disconnect()
    else:
        print("\n❌ Failed to connect")
