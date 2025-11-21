/**
 * Progress Bar Component - Top loading bar for page transitions
 */
import React, { useEffect, useState } from 'react'
import { useLocation } from 'react-router-dom'

export const ProgressBar: React.FC = () => {
  const location = useLocation()
  const [progress, setProgress] = useState(0)
  const [isLoading, setIsLoading] = useState(false)

  useEffect(() => {
    setIsLoading(true)
    setProgress(20)

    const timer1 = setTimeout(() => setProgress(40), 100)
    const timer2 = setTimeout(() => setProgress(60), 200)
    const timer3 = setTimeout(() => setProgress(80), 300)
    const timer4 = setTimeout(() => {
      setProgress(100)
      setTimeout(() => {
        setIsLoading(false)
        setProgress(0)
      }, 200)
    }, 400)

    return () => {
      clearTimeout(timer1)
      clearTimeout(timer2)
      clearTimeout(timer3)
      clearTimeout(timer4)
    }
  }, [location])

  if (!isLoading && progress === 0) return null

  return (
    <div
      style={{
        position: 'fixed',
        top: 0,
        left: 0,
        right: 0,
        height: 3,
        zIndex: 9999,
        background: 'linear-gradient(90deg, #2563eb, #1d4ed8)',
        width: `${progress}%`,
        transition: 'width 0.2s ease-out, opacity 0.2s ease-out',
        opacity: progress === 100 ? 0 : 1,
        boxShadow: '0 0 10px rgba(37, 99, 235, 0.5)',
      }}
    />
  )
}

export default ProgressBar
