import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { nodePolyfills } from 'vite-plugin-node-polyfills'
import { fileURLToPath, URL } from 'node:url'

// Builds Ketcher 3.18.0 into the Django static tree. The bundle is emitted as
// a single ES module and loaded with <script type="module">.
//
// Output:
//   static/js/ketcher3/ketcher.bundle.js   (ESM: React + ketcher-react + ketcher-standalone)
//   static/js/ketcher3/styles.css          (ketcher-react styles, extracted by Vite lib mode)
//   static/js/ketcher3/assets/...          (fonts/images/wasm emitted by Ketcher, if any)
export default defineConfig({
  plugins: [
    react(),
    // Ketcher (indigo) imports Node built-ins such as `events` (EventEmitter).
    // Vite does not polyfill Node built-ins by default; this plugin provides
    // browser polyfills for events/buffer/process/util/stream/etc.
    nodePolyfills(),
  ],
  build: {
    outDir: fileURLToPath(new URL('../../static/js/ketcher3', import.meta.url)),
    emptyOutDir: true,
    // ketcher-core ships ESM chunks that still contain conditional
    // `require("raphael")` calls (guarded by `typeof window`). By default the
    // CommonJS plugin skips files that mix import/export with require, leaving
    // a bare `require` in the browser bundle. This forces those to be converted.
    commonjsOptions: {
      transformMixedEsModules: true,
    },
    lib: {
      entry: fileURLToPath(new URL('src/main.jsx', import.meta.url)),
      name: 'KetcherApp',
      // ESM (not IIFE) so Rollup does not emit bare `require(...)` calls from
      // CommonJS dependencies (which break in the browser: "require is not defined").
      formats: ['es'],
      fileName: () => 'ketcher.bundle.js',
    },
    rollupOptions: {
      output: {
        // Ketcher uses dynamic imports internally; inlining them keeps the
        // whole editor in a single file.
        inlineDynamicImports: true,
      },
    },
  },
})
