import path from "path";
import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";
import { Buffer } from 'buffer';

export default defineConfig({
  plugins: [react()],
  logLevel: 'warn', // Reduce build output noise (suppresses 'use client' warnings)
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
      'buffer': 'buffer',
      'path': 'path-browserify',
      '~': path.resolve(__dirname, "./src"),
      '#minpath': path.resolve(__dirname, "./node_modules/vfile/lib/minpath.js"),
      '#minproc': path.resolve(__dirname, "./node_modules/vfile/lib/minproc.js"),
      '#minurl': path.resolve(__dirname, "./node_modules/vfile/lib/minurl.js"),
      'node:url': 'url',
    },
  },
  define: {
    'global.Buffer': Buffer,
    global: "globalThis",
    "process.env": {},
  },
  base: "/",
  server: {
    // No proxy needed for standalone frontend deployment
  },
  build: {
    // Increase chunk size warning limit
    chunkSizeWarningLimit: 1000,
    rollupOptions: {
      external: [],
      output: {
        // Manual chunk splitting to reduce memory usage
        manualChunks: (id) => {
          // Split node_modules into separate chunks to reduce memory usage
          if (id.includes('node_modules')) {
            // Large libraries get their own chunks
            if (id.includes('@mui')) {
              return 'mui';
            }
            if (id.includes('react') || id.includes('react-dom')) {
              return 'react-vendor';
            }
            if (id.includes('@supabase')) {
              return 'supabase';
            }
            if (id.includes('recharts')) {
              return 'recharts';
            }
            if (id.includes('@google/generative-ai')) {
              return 'google-ai';
            }
            if (id.includes('openai')) {
              return 'openai';
            }
            if (id.includes('ethers') || id.includes('@thirdweb')) {
              return 'web3';
            }
            // Split remaining vendor into smaller chunks
            if (id.includes('node_modules')) {
              // Group by first letter to create smaller chunks
              const match = id.match(/node_modules\/(@[^/]+\/)?([^/]+)/);
              if (match) {
                const pkgName = match[2] || match[1] || '';
                if (pkgName.startsWith('@')) {
                  return 'vendor-scoped';
                }
                const firstLetter = pkgName.charAt(0).toLowerCase();
                if (firstLetter >= 'a' && firstLetter <= 'm') {
                  return 'vendor-a-m';
                }
                return 'vendor-n-z';
              }
            }
            return 'vendor';
          }
        },
      },
    },
    // Reduce memory usage during build
    minify: 'esbuild', // Faster than terser
    target: 'es2015', // Broader compatibility, less transformation
  },
});
