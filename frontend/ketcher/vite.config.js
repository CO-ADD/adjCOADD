import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { fileURLToPath, URL } from 'node:url'

// Builds Ketcher 3.18.0 into the Django static tree so it can be loaded with a
// plain <script src> / <link> (the same way the old static/js/ketcher bundle was).
//
// Output:
//   static/js/ketcher3/ketcher.bundle.js   (IIFE: React + ketcher-react + ketcher-standalone)
//   static/js/ketcher3/ketcher.bundle.css  (ketcher-react styles, extracted)
//   static/js/ketcher3/assets/...          (fonts/images/wasm emitted by Ketcher, if any)
export default defineConfig({
  plugins: [react()],
  build: {
    outDir: fileURLToPath(new URL('../../static/js/ketcher3', import.meta.url)),
    emptyOutDir: true,
    lib: {
      entry: fileURLToPath(new URL('src/main.jsx', import.meta.url)),
      name: 'KetcherApp',
      formats: ['iife'],
      fileName: () => 'ketcher.bundle.js',
    },
    rollupOptions: {
      output: {
        // Ketcher uses dynamic imports internally; a single-file IIFE bundle
        // requires inlining them.
        inlineDynamicImports: true,
      },
    },
  },
})
