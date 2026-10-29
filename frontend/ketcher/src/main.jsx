import React, { useState } from 'react'
import { createRoot } from 'react-dom/client'
import { StandaloneStructServiceProvider } from 'ketcher-standalone'
import { Editor } from 'ketcher-react'
import 'ketcher-react/dist/index.css'

// Standalone mode: the indigo struct service is bundled with Ketcher, so no
// remote struct-service URL is required.
const structServiceProvider = new StandaloneStructServiceProvider()

function KetcherEditor({ initialMol }) {
  const [ketcher, setKetcher] = useState(null)

  function handleOnInit(k) {
    setKetcher(k)
    // Expose the instance so page scripts can call getMolfile()/setMolecule().
    window.ketcher = k
    if (initialMol) {
      k.setMolecule(initialMol.replace(/\\n/g, '\n'))
    }
  }

  function updateSmiles() {
    ketcher.getMolfile('v2000').then((molfile) => {
      const out = document.getElementById('id_smol')
      if (out) out.value = molfile
    })
  }

  return (
    <div>
      <div>
        <button id="get_new_mol" type="button" onClick={updateSmiles}>
          Update Smiles
        </button>
      </div>
      <Editor
        staticResourcesUrl=""
        structServiceProvider={structServiceProvider}
        onInit={handleOnInit}
      />
    </div>
  )
}

// Same contract as the old static/js/ketcher bundle:
//   - a hidden <input id="mol"> carries the initial molfile
//   - <div id="root"> is the mount point
//   - the current molfile is written to <input id="id_smol"> on "Update Smiles"
//   - window.ketcher is exposed for external Check/Save handlers
const rootEl = document.getElementById('root')
if (rootEl) {
  const molInput = document.getElementById('mol')
  const initialMol = molInput ? molInput.value : ''
  createRoot(rootEl).render(<KetcherEditor initialMol={initialMol} />)
}
